#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

set -e
set -u
set -o pipefail

####################################################################
# ninjaMap.sh
#
# Accurate alignment of reads from synthetic microbiome to a comebined
# set of microbial genomes using bowtie2. The accuracy comes from taking
# care of reads that map equally well to multiple references.
#
# This script is used in conjunction with aegea batch
#
# Variables required from sourcing script
# coreN=4; numPerCore=1G; maxInsert=3000; maxAlignments=200;
# S3DBPATH=/czbiohub-microbiome/Synthetic_Community/Genome_References/Bowtie2Index_090718
# SAMPLE_NAME; fastq1=/czbiohub-microbiome/...; fastq2;
# BAM_OUTPUT; REL_AB_OUTPUT; READ_ACC_OUTPUT
#
# bbtools assumes 16G of memory -Xmx16g; needs sambamba in conda env
######################################################################
START_TIME=$SECONDS
# export $@
export PATH="/opt/conda/bin:${PATH}"

sampleRate="${sampleRate:-1}"
coreNum="${coreNum:-16}"
coreN="${coreN:-15}"
memPerCore="${memPerCore:-2G}"
maxInsert="${maxInsert:-3000}"
maxAlignments="${maxAlignments:-200}"
minPercId="${minPercId:-0}"
minReadQuality="${minReadQuality:-0}"
minMapQuality="${minMapQuality:-10}"
minAlnCov="${minAlnCov:-0}"
export NUMEXPR_MAX_THREADS="${coreN}" # required for numpy
human="${human:-1}"
mouse="${mouse:-1}"
# Inputs
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/2019-05-16_StrainVerification/Dorea-longicatena-DSM-13814
# fastq1=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz

# S3DBPATH="s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/20180725/scv1/db/"
# REFDBNAME="uniform10x_ninjaIndex.ninjaIndex"
# BINMAP_FILENAME="uniform10x_ninjaIndex.ninjaIndex.binmap.csv"

S3DBPATH=${S3DBPATH:-"s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Narrow/20190911/scv2/db/"}
REFDBNAME=${REFDBNAME:-"20190911_scv2"}
STRAIN_MAP_FILENAME=${STRAIN_MAP_FILENAME:-"20190911_scv2.ninjaIndex.binmap.csv"}

SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

echo $PATH
LOCAL=$(pwd)

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
DOWNLOAD_DB="${OUTPUTDIR}/download_db"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
TMP_OUTPUTS="${OUTPUTDIR}/bowtie2"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
STATS_DIR="${LOCAL_OUTPUT}/Stats"
BOWTIE2_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
NINJA_OUTPUT="${LOCAL_OUTPUT}/ninjaMap"
GENOME_COV_OUTPUT="${LOCAL_OUTPUT}/genome_coverage"
S3OUTPUTPATH=${S3OUTPUTPATH%/}
S3DBPATH=${S3DBPATH%/}
LOCAL_DB_PATH="${OUTPUTDIR}/reference"
OUTPUT_PREFIX="${SAMPLE_NAME}.sortedByCoord"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}" "${RAW_FASTQ}" "${QC_FASTQ}" "${STATS_DIR}" "${DOWNLOAD_DB}"
mkdir -p "${LOCAL_DB_PATH}" "${BOWTIE2_OUTPUT}" "${NINJA_OUTPUT}" "${TMP_OUTPUTS}" "${GENOME_COV_OUTPUT}"

trap '{aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}";
    rm -rf ${OUTPUTDIR} ;
    exit 255; }' 1

adapterFile="adapters,phix"
# offLimitRegions="./data/combined_excluded_regions_threshold9.bed"
scriptFolder="./scripts"
BOWTIE2_DB=${LOCAL_DB_PATH}/bowtie2_index/${REFDBNAME}
# REF_FASTA=${LOCAL_DB_PATH}/${REFDBNAME}.fna

# Copy genome reference over
aws s3 sync --quiet ${S3DBPATH}/ ${LOCAL_DB_PATH}/
referenceNameFile=${LOCAL_DB_PATH}/${STRAIN_MAP_FILENAME}

# Constant definitions for bbduk
trimQuality="${trimQuality:-25}"
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}

echo "Starting to Process Sample: "${SAMPLE_NAME}

# Copy fastq.gz files from S3, only 2 files per sample
aws s3 cp --quiet ${fastq1} ${RAW_FASTQ}/read1.fastq.gz
aws s3 cp --quiet ${fastq2} ${RAW_FASTQ}/read2.fastq.gz

# Downsample reads to get results faster
reformat.sh \
samplerate=${sampleRate} \
sampleseed=1772 \
in="${RAW_FASTQ}/read1.fastq.gz" in2="${RAW_FASTQ}/read2.fastq.gz" \
out="${RAW_FASTQ}/read1_sampled.fastq.gz" out2="${RAW_FASTQ}/read2_sampled.fastq.gz" \
&> ${LOG_DIR}/reformat.log.txt

# Use bbduk to trim reads, -eoom exits when out of memory
bbduk.sh -Xmx60g -eoom \
    in1=${RAW_FASTQ}/read1_sampled.fastq.gz \
    in2=${RAW_FASTQ}/read2_sampled.fastq.gz \
    out1=${QC_FASTQ}/read1_trimmed.fastq.gz \
    out2=${QC_FASTQ}/read2_trimmed.fastq.gz \
    ref=${adapterFile} \
    ktrim=r \
    k=${kmer_value} \
    mink=${min_kmer_value} \
    hdist=1 tbo qtrim=rl \
    trimq=${trimQuality} \
    minlen=${minLength} \
    refstats=${STATS_DIR}/adapter_trimming_stats_per_ref.txt \
    &>> ${LOG_DIR}/bbduk.log.txt

# ref: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
if [[ ${human} -eq 1 ]]
then
    # download genome index
    aws s3 sync --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/human/ ${DOWNLOAD_DB}/GRCh38/
    #aws s3 cp --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/GRCh38_noalt_as.zip ${DOWNLOAD_DB}/GRCh38_noalt_as.zip
    #unzip "${DOWNLOAD_DB}/GRCh38_noalt_as.zip"
    # bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
    bowtie2 -p ${coreN} -x "${DOWNLOAD_DB}/GRCh38/human" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" -S "${QC_FASTQ}/Human_mapped_and_unmapped.sam"
    #convert file .sam to .bam
    samtools view -bS "${QC_FASTQ}/Human_mapped_and_unmapped.sam" > "${QC_FASTQ}/Human_mapped_and_unmapped.bam"
    samtools flagstat "${QC_FASTQ}/Human_mapped_and_unmapped.bam" > "${LOG_DIR}/Human_contamination-stat.txt"
    # SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
    samtools view -b -f 12 -F 256 "${QC_FASTQ}/Human_mapped_and_unmapped.bam" > "${QC_FASTQ}/Human_bothReadsUnmapped.bam"
    # sort bam file by read name (-n) to have paired reads next to each other
    samtools sort -n -m ${memPerCore} -@ ${coreN} "${QC_FASTQ}/Human_bothReadsUnmapped.bam" -o "${QC_FASTQ}/Human_bothReadsUnmapped_sorted.bam"
    samtools fastq -@ ${coreN} "${QC_FASTQ}/Human_bothReadsUnmapped_sorted.bam" \
    -1 "${BOWTIE2_OUTPUT}/Human_host_removed.read1.fastq.gz" \
    -2 "${BOWTIE2_OUTPUT}/Human_host_removed.read2.fastq.gz" \
    -0 /dev/null -s /dev/null -n
  fi

if [[ ${mouse} -eq 1 ]]
then
  # download genome index
  aws s3 sync --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/mouse/ ${DOWNLOAD_DB}/GRCm39/
  #aws s3 cp --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/GRCm39.zip ${DOWNLOAD_DB}/GRCm39.zip
  #unzip "${DOWNLOAD_DB}/GRCm39.zip"
  # bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
  bowtie2 -p ${coreN} -x "${DOWNLOAD_DB}/GRCm39/mouse" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" -S "${QC_FASTQ}/Mouse_mapped_and_unmapped.sam"
  #convert file .sam to .bam
  samtools view -bS "${QC_FASTQ}/Mouse_mapped_and_unmapped.sam" > "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam"
  samtools flagstat "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam" > "${LOG_DIR}/Mouse_contamination-stat.txt"
  # SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
  samtools view -b -f 12 -F 256 "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam" > "${QC_FASTQ}/Mouse_bothReadsUnmapped.bam"
  # sort bam file by read name (-n) to have paired reads next to each other
  samtools sort -n -m ${memPerCore} -@ ${coreN} "${QC_FASTQ}/Mouse_bothReadsUnmapped.bam" -o "${QC_FASTQ}/Mouse_bothReadsUnmapped_sorted.bam"
  samtools fastq -@ ${coreN} "${QC_FASTQ}/Mouse_bothReadsUnmapped_sorted.bam" \
  -1 "${BOWTIE2_OUTPUT}/Mouse_host_removed.read1.fastq.gz" \
  -2 "${BOWTIE2_OUTPUT}/Mouse_host_removed.read2.fastq.gz" \
  -0 /dev/null -s /dev/null -n
fi

echo "SAMPLE_NAME,Host_Contaminant,Total_mapped_reads,Mapped_rate(%),Total_mapped_paired_reads,Mapped_Paired_rate(%)" >> ${LOG_DIR}/Host_Contaminants_stats.csv

LINE="${SAMPLE_NAME},Human,"
mapped_reads=`grep "mapped ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d " " -f 1`
mapped_rate=`grep "mapped ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
mapped_paired_reads=`grep "paired ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d " " -f 1`
mapped_paired_rate=`grep "paired ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
LINE+="${mapped_reads},${mapped_rate},${mapped_paired_reads},${mapped_paired_rate}"
echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv

LINE="${SAMPLE_NAME},Mouse,"
mapped_reads=`grep "mapped ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d " " -f 1`
mapped_rate=`grep "mapped ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
mapped_paired_reads=`grep "paired ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d " " -f 1`
mapped_paired_rate=`grep "paired ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
LINE+="${mapped_reads},${mapped_rate},${mapped_paired_reads},${mapped_paired_rate}"
echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv


:<< 'COMMENT'
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
bowtie2 \
    --very-sensitive \
    -X ${maxInsert} \
    -k ${maxAlignments} \
    --threads ${coreN} \
    -x ${BOWTIE2_DB} \
    --no-mixed \
    --no-discordant \
    --end-to-end \
    --no-unal \
    -1 ${QC_FASTQ}/read1_trimmed.fastq.gz \
    -2 ${QC_FASTQ}/read2_trimmed.fastq.gz | \
    samtools view \
        -@ ${coreN} \
        -bh \
        -o ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam - |\
    tee -a ${LOG_DIR}/read_mapping.log.txt

# Original bowtie2 parameters
# Removed: -f 3 \
# Removed: -D 10 -R 2 -L 31 -i S,0,2.50 -N 0
# Added: --very-sensitive
# 20200504: (based on: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#mates-can-overlap-contain-or-dovetail-each-other)
#           Try adding --no-overlap and --no-contain. Since this should reduce the number of spurious matches, also try
#           replacing -k ${maxAlignments} with -a for all.

# Fix Mates
samtools sort \
    -n \
    -@ ${coreN} \
    -m ${memPerCore} \
    -T ${OUTPUTDIR} \
    -O BAM \
    ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam |\
samtools fixmate \
    -O BAM \
    -cm \
    - - | \
samtools sort \
    -@ ${coreN} \
    -m ${memPerCore} \
    -T ${OUTPUTDIR} \
    -o ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam

# samtools sort -@ ${coreN} -o ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam
samtools index -@ ${coreN} ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam

# rmove unsorted bam file to save space
rm ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam
# 3328 =
#   not primary alignment (0x100)
#   read is PCR or optical duplicate (0x400)
#   supplementary alignment (0x800)
# samtools view -F 3328 -q 10 Dorea-longicatena-DSM-13814.processed.bam | cut -f1 | sort | uniq | wc -l
#  -threads ${coreN} \     -threads ${coreNum} \
timem python ${scriptFolder}/ninjaMap_parallel_2.py \
    -bam ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam \
    -bin ${referenceNameFile} \
    -outdir ${NINJA_OUTPUT} \
    -prefix ${SAMPLE_NAME} |\
    tee -a ${LOG_DIR}/${SAMPLE_NAME}.ninjaMap.log.txt

# Tabulate read count
totalReads=$(( $( zcat ${RAW_FASTQ}/read1_sampled.fastq.gz | wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat ${QC_FASTQ}/read1_trimmed.fastq.gz | wc -l ) / 4 ))
uniqueReads=$( samtools view -f 0x40 ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam | cut -f1 | sort -u | wc -l )
echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned' > ${STATS_DIR}/read_accounting.csv
echo ${SAMPLE_NAME}','${totalReads}','${readsAfterTrim}','${uniqueReads} >> ${STATS_DIR}/read_accounting.csv

COMMENT

echo "Contamination check completed."
ls ${LOCAL}
du -sh ${LOCAL}
date
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job_contamination.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job_contamination.complete
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
