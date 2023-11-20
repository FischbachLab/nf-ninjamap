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
coreN=$coreNum
memPerCore="${memPerCore:-2G}"
maxInsert="${maxInsert:-3000}"
maxAlignments="${maxAlignments:-200}"
minPercId="${minPercId:-0}"
minReadQuality="${minReadQuality:-0}"
minMapQuality="${minMapQuality:-10}"
minAlnCov="${minAlnCov:-0}"
export NUMEXPR_MAX_THREADS="${coreN}" # required for numpy
human="${human:-0}"
mouse="${mouse:-0}"
coverage="${coverage:-0}"
# Inputs
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/2019-05-16_StrainVerification/Dorea-longicatena-DSM-13814
# fastq1=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R1_001.fastq.gz
# fastq2=s3://czbiohub-microbiome/Original_Sequencing_Data/180727_A00111_0179_BH72VVDSXX/Alice_Cheng/Strain_Verification/Dorea-longicatena-DSM-13814_S275_R2_001.fastq.gz

# S3DBPATH="s3://czbiohub-microbiome/ReferenceDBs/NinjaMap/Index/20180725/scv1/db/"
# REFDBNAME="uniform10x_ninjaIndex.ninjaIndex"
# BINMAP_FILENAME="uniform10x_ninjaIndex.ninjaIndex.binmap.csv"

S3DBPATH=${S3DBPATH:-"s3://maf-versioned/ninjamap/Index/HCom2/db/"}
REFDBNAME=${REFDBNAME:-"HCom2"}
STRAIN_MAP_FILENAME=${STRAIN_MAP_FILENAME:-"HCom2.ninjaIndex.binmap.csv"}

SAMPLE_NAME=$(basename ${S3OUTPUTPATH})

echo $PATH
LOCAL=.

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )
RAW_FASTQ="${OUTPUTDIR}/raw_fastq"
DOWNLOAD_DB="${OUTPUTDIR}/download_db"
QC_FASTQ="${OUTPUTDIR}/trimmed_fastq"
TMP_OUTPUTS="${OUTPUTDIR}/bowtie2"
GENOME_COV_OUTPUT="${OUTPUTDIR}/genome_coverage"
LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
STATS_DIR="${LOCAL_OUTPUT}/Stats"
BOWTIE2_OUTPUT="${LOCAL_OUTPUT}/bowtie2"
NINJA_OUTPUT="${LOCAL_OUTPUT}/ninjaMap"
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
scriptFolder="/work/scripts"
BOWTIE2_DB=${LOCAL_DB_PATH}/bowtie2_index/${REFDBNAME}
# REF_FASTA=${LOCAL_DB_PATH}/${REFDBNAME}.fna

# Copy genome reference over  ${params.db_path}/${params.db}/db/
# Check the reference location
if [[ $S3DBPATH = s3* ]]
then
	echo "Sync index files from aws s3"
  aws s3 sync --quiet ${S3DBPATH}/ ${LOCAL_DB_PATH}/
else
  echo "Copying index files from the external URL"
  #S3DBPATH=${S3DBPATH%%/HCom2*}
  S3DBPATH="https://zenodo.org/record/7872423/files/hCom2_20221117.ninjaIndex.tar.gz"
  wget $S3DBPATH -P ${LOCAL_DB_PATH}/
  tar -xzvf ${LOCAL_DB_PATH}/hCom2_20221117.ninjaIndex.tar.gz  -C ${LOCAL_DB_PATH}/
  mv ${LOCAL_DB_PATH}/HCom2_20221117/db/* ${LOCAL_DB_PATH}/
fi

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

: <<'END'
bbmap not working
File "/mnt/efs/scratch/Xmeng/data/ninjaMap/ninjaMap_555.py", line 556, in is_perfect_alignment
  edit_dist = dict(aln.tags)['NM']
KeyError: 'NM'

    bbmap.sh \
    perfectmode=t ambiguous=all  pairlen=3000 \
    nmtag=t \
    in=${RAW_FASTQ}/read1_sampled.fastq.gz  \
    in2=${RAW_FASTQ}/read2_sampled.fastq.gz  \
    ref="${LOCAL_DB_PATH}/HCom2.fna" \
    out=stdout | samtools view -@ ${coreN} -bh -o ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam -

cp ${RAW_FASTQ}/read1_sampled.fastq.gz ${QC_FASTQ}/read1_trimmed.fastq.gz
cp ${RAW_FASTQ}/read2_sampled.fastq.gz ${QC_FASTQ}/read2_trimmed.fastq.gz
END

# Try to use all available cores to speed up alignments
#mycpu=`grep -c ^processor /proc/cpuinfo`
mycpu=$(bc <<< "(`grep -c ^processor /proc/cpuinfo` - 1)" )
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
bowtie2 \
    --very-sensitive \
    -X ${maxInsert} \
    -k ${maxAlignments} \
    --threads ${mycpu} \
    -x ${BOWTIE2_DB} \
    --no-mixed \
    --no-discordant \
    --end-to-end \
    --no-unal \
    -1 ${QC_FASTQ}/read1_trimmed.fastq.gz \
    -2 ${QC_FASTQ}/read2_trimmed.fastq.gz | \
    samtools view \
        -@ ${mycpu} \
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
 #-m ${memPerCore} \   #-m ${memPerCore} \
# Fix Mates
samtools sort \
    -n \
    -@ ${coreN} \
    -T ${OUTPUTDIR} \
    -O BAM \
    ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam |\
samtools fixmate \
    -O BAM \
    -cm \
    - - | \
samtools sort \
    -@ ${coreN} \
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
# check if adding -coverage option
if [ ${coverage} -eq 1 ];
then
  cov_option="-coverage"
else
  cov_option=""
fi

# get bam file size
bamSizeCheck=$(stat -c %s "${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam")

all_mapped_reads=2
# 100k?
if [ $bamSizeCheck -lt 100000 ]; then
    samtools flagstat "${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam" > "${LOG_DIR}/Alignment-stat.txt"
    all_mapped_reads=`grep "mapped ("  "${LOG_DIR}/Alignment-stat.txt" | cut -d " " -f 1`
fi

# make sure the number of aligned reads > 1
if [ $all_mapped_reads -gt 1 ]; then
python ${scriptFolder}/ninjaMap_parallel_5.py \
    -bam ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam \
    -bin ${referenceNameFile} \
    -outdir ${NINJA_OUTPUT} \
    ${cov_option} \
    -prefix ${SAMPLE_NAME} |\
    tee -a ${LOG_DIR}/${SAMPLE_NAME}.ninjaMap.log.txt
else
    rm -rf ${OUTPUTDIR}
    exit 0
fi

# Tabulate read count
totalReads=$(( $( zcat ${RAW_FASTQ}/read1_sampled.fastq.gz | wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat ${QC_FASTQ}/read1_trimmed.fastq.gz | wc -l ) / 4 ))
uniqueReads=$( samtools view -f 0x40 ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam | cut -f1 | sort -u | wc -l )
echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned' > ${STATS_DIR}/read_accounting.csv
echo ${SAMPLE_NAME}','${totalReads}','${readsAfterTrim}','${uniqueReads} >> ${STATS_DIR}/read_accounting.csv

#################################################################
#Add the overall Percent_Coverage to the abundance table
#################################################################
#filter bam file for perfect alignment
if [ $all_mapped_reads -gt 1 ]; then
python ${scriptFolder}/bamParser.py -bam ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam -id 100 -aln_len 100 -out ${GENOME_COV_OUTPUT}/${SAMPLE_NAME}.filtered.bam
pileup.sh in=${GENOME_COV_OUTPUT}/${SAMPLE_NAME}.filtered.bam out=${GENOME_COV_OUTPUT}/${SAMPLE_NAME}_coverage.txt overwrite=t delcoverage=f 2>${GENOME_COV_OUTPUT}/${SAMPLE_NAME}_stats.txt

samtools faidx ${LOCAL_DB_PATH}/${REFDBNAME}.fna && cut -f1,2 ${LOCAL_DB_PATH}/${REFDBNAME}.fna.fai > ${GENOME_COV_OUTPUT}/Ninja.genomes
cut -f1  ${GENOME_COV_OUTPUT}/Ninja.genomes | sed  's/_Node.*//' | sort | uniq  > ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt
echo -e "Stain_Name\tPercent_Coverage\tCoverage_Depth" > ${NINJA_OUTPUT}/${SAMPLE_NAME}_summary_coverage.tsv
for i in $(cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt)
do
  grep $i ${GENOME_COV_OUTPUT}/${SAMPLE_NAME}_coverage.txt | awk '{c+=$6; s+=$2*$3;l+=$3}END{print $1"\t"c/l*100"\t"s/l}' >> ${NINJA_OUTPUT}/${SAMPLE_NAME}_summary_coverage.tsv
done
cut -f2,3 ${NINJA_OUTPUT}/${SAMPLE_NAME}_summary_coverage.tsv | paste ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3}' > ${NINJA_OUTPUT}/tmp0.ninjaMap.abundance.csv

mv ${NINJA_OUTPUT}/tmp0.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv
rm ${NINJA_OUTPUT}/${SAMPLE_NAME}_summary_coverage.tsv

fi
#################################################################
#Add Singular and Escrow coverage depth to the abundance table
#Notet that the sum of Singular and Escrow coverage coverage is often over 100%
################################################################

if [ ${coverage} -eq 1 ];
then
    s="singular"
    samtools sort -@ ${coreN} ${NINJA_OUTPUT}/${SAMPLE_NAME}.singular.bam > ${NINJA_OUTPUT}/${SAMPLE_NAME}.singular_sorted.bam
    bedtools bamtobed -i ${NINJA_OUTPUT}/${SAMPLE_NAME}.singular_sorted.bam > ${NINJA_OUTPUT}/${s}.bed
    echo -e "Strain_Name\tSingular_Coverage\tSingular_Depth\tSingular_Bases" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv
    #for i in $(cut -f1 Ninja.genomes | sed  's/_Node.*//' |  uniq )  | cut -f1,2,3 -
    for i in $(cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt)
    do
        if [ $(grep -c "^$i" ${NINJA_OUTPUT}/${s}.bed) -eq 0 ];
        then
          > ${GENOME_COV_OUTPUT}/tmp.bed
        else
          grep "^$i" ${NINJA_OUTPUT}/${s}.bed | sort | uniq | cut -f1,2,3  > ${GENOME_COV_OUTPUT}/tmp.bed
        fi
        grep "^$i"  ${GENOME_COV_OUTPUT}/Ninja.genomes > ${GENOME_COV_OUTPUT}/tmp.genome
        bedtools genomecov -i ${GENOME_COV_OUTPUT}/tmp.bed -g ${GENOME_COV_OUTPUT}/tmp.genome -bga > ${GENOME_COV_OUTPUT}/tmps.bedg
        #compute coverage
        zero=$(awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
        nonzero=$( awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
        if [ -z $zero ];  then zero=0; fi
        if [ -z $nonzero ];  then nonzero=0; fi
        cov_pct=$(bc <<< "scale=4; ($nonzero / ($zero + $nonzero))*100")
        # compute depth
        totalbases=$(awk '{t+=$4*($3-$2)} END {print t}' ${GENOME_COV_OUTPUT}/tmps.bedg)
        if [ $nonzero -eq 0 ] || [ $totalbases -eq 0 ] ;
        then
           cov_depth=0
        else
           cov_depth=$(bc <<< "scale=4; ($totalbases / $nonzero)")
        fi
        printf "${i}\t${cov_pct}\t${cov_depth}\t${nonzero}\n" >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
     done
     cut -f2,3,4 ${NINJA_OUTPUT}/${s}_summary_depth.tsv | paste ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3","$4}' > ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv

     s="escrow"
     samtools sort -@ ${coreN} ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow.bam > ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam
     bedtools bamtobed -i ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam  > ${NINJA_OUTPUT}/${s}.bed
     echo -e "Strain_Name\tEscrow_Coverage\tEscorw_Depth\tEscorw_Bases" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv
     for i in $(cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt)
     do
         if [ $(grep -c "^$i" ${NINJA_OUTPUT}/${s}.bed) -eq 0 ];
         then
            > ${GENOME_COV_OUTPUT}/tmp.bed
         else
           grep "^$i" ${NINJA_OUTPUT}/${s}.bed | sort | uniq | cut -f1,2,3  > ${GENOME_COV_OUTPUT}/tmp.bed
         fi
         grep "^$i"  ${GENOME_COV_OUTPUT}/Ninja.genomes > ${GENOME_COV_OUTPUT}/tmp.genome
         bedtools genomecov -i ${GENOME_COV_OUTPUT}/tmp.bed -g ${GENOME_COV_OUTPUT}/tmp.genome -bga > ${GENOME_COV_OUTPUT}/tmps.bedg
         #compute coverage
         zero=$(awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         nonzero=$( awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         if [ -z $zero ];  then zero=0; fi
         if [ -z $nonzero ];  then nonzero=0; fi
         cov_pct=$(bc <<< "scale=4; ($nonzero / ($zero + $nonzero))*100")
         # compute depth
         totalbases=$(awk '{t+=$4*($3-$2)} END {print t}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         if  [ $nonzero -eq 0 ] || [ $totalbases -eq 0 ];
         then
            cov_depth=0
         else
            cov_depth=$(bc <<< "scale=4; ($totalbases / $nonzero)")
         fi
         printf "${i}\t${cov_pct}\t${cov_depth}\t${nonzero}\n" >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
      done
      cut -f2,3,4 ${NINJA_OUTPUT}/${s}_summary_depth.tsv | paste ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3","$4}' > ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv
      mv ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv
      rm ${NINJA_OUTPUT}/*_summary_depth.tsv ${NINJA_OUTPUT}/tmp*.ninjaMap.abundance.csv
      rm ${NINJA_OUTPUT}/${SAMPLE_NAME}.singular.bam  ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow.bam
      rm ${NINJA_OUTPUT}/*.bed
fi

:<<"COMM"
if [ ${coverage} -eq 1 ];
then
    s="singular"
    pileup.sh in=${NINJA_OUTPUT}/${SAMPLE_NAME}.singular.bam out=${NINJA_OUTPUT}/${s}_coverage.txt overwrite=t delcoverage=f 2>${NINJA_OUTPUT}/${s}_stats.txt
    echo -e "Stain_Name\tSingular_Depth\tSingular_Coverage" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv
    #cut -f1 ${NINJA_OUTPUT}/${s}_coverage.txt | tail -n+2  | sed  's/_Node.*//'
    for i in $(cut -f1 ${NINJA_OUTPUT}/${s}_coverage.txt | tail -n+2  | sed  's/_Node.*//' |  uniq )
    do
     grep $i ${NINJA_OUTPUT}/${s}_coverage.txt | awk '{t+=$2*$3; c+=$6; l+=$3}END{print $1"\t"t/l"\t"c/l}' >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
    done
    cut -f2,3 ${NINJA_OUTPUT}/${s}_summary_depth.tsv | paste ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3}' > ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv

    s="escrow"
    pileup.sh in=${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow.bam out=${NINJA_OUTPUT}/${s}_coverage.txt overwrite=t delcoverage=f 2>${NINJA_OUTPUT}/${s}_stats.txt
    echo -e "Stain_Name\tEscrow_Depth\tEscorw_Coverage" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv
    for i in $(cut -f1 ${NINJA_OUTPUT}/${s}_coverage.txt | tail -n+2  | sed  's/_Node.*//' |  uniq )
    do
     grep $i ${NINJA_OUTPUT}/${s}_coverage.txt | awk '{t+=$2*$3; c+=$6; l+=$3}END{print $1"\t"t/l"\t"c/l}' >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
    done
    cut -f2,3 ${NINJA_OUTPUT}/${s}_summary_depth.tsv | paste ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3}' > ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv

    mv ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv
    rm ${NINJA_OUTPUT}/*_summary_depth.tsv  ${NINJA_OUTPUT}/*_coverage.txt ${NINJA_OUTPUT}/*_stats.txt  ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv
fi
COMM
######################################################################
#  HOST contamination check
######################################################################
alignmentRate=$(( uniqueReads * 100 / readsAfterTrim ))
#if [ ${alignmentRate} -lt 95 ] || [ ${human} -eq 1 ] || [ ${mouse} -eq 1 ];
#then
echo "SAMPLE_NAME,Host_Contaminant,Total_mapped_reads,Mapped_rate(%),Total_mapped_paired_reads,Mapped_Paired_rate(%)" >> ${LOG_DIR}/Host_Contaminants_stats.csv
#fi

# ref: https://www.metagenomics.wiki/tools/short-read/remove-host-sequences
if [ ${alignmentRate} -lt 95 ] || [ ${human} -eq 1 ];
then
    # download genome index
    #aws s3 sync --quiet s3://maf-versioned/ninjamap/bowtie2/human/ ${DOWNLOAD_DB}/GRCh38/
    # link genome index
    ln -s /mnt/efs/databases/Bowtie2/human ${DOWNLOAD_DB}/GRCh38
    #aws s3 cp --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/GRCh38_noalt_as.zip ${DOWNLOAD_DB}/GRCh38_noalt_as.zip
    #unzip "${DOWNLOAD_DB}/GRCh38_noalt_as.zip"
    # bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
    bowtie2 -p ${mycpu} -x "${DOWNLOAD_DB}/GRCh38/human" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" -S "${QC_FASTQ}/Human_mapped_and_unmapped.sam"
    #convert file .sam to .bam
    samtools view -bS "${QC_FASTQ}/Human_mapped_and_unmapped.sam" > "${QC_FASTQ}/Human_mapped_and_unmapped.bam"
    samtools flagstat "${QC_FASTQ}/Human_mapped_and_unmapped.bam" > "${LOG_DIR}/Human_contamination-stat.txt"
    # SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
    samtools view -b -f 12 -F 256 "${QC_FASTQ}/Human_mapped_and_unmapped.bam" > "${QC_FASTQ}/Human_bothReadsUnmapped.bam"
    # sort bam file by read name (-n) to have paired reads next to each other
    samtools sort -n -@ ${mycpu} "${QC_FASTQ}/Human_bothReadsUnmapped.bam" -o "${QC_FASTQ}/Human_bothReadsUnmapped_sorted.bam"
    samtools fastq -@ ${mycpu} "${QC_FASTQ}/Human_bothReadsUnmapped_sorted.bam" \
    -1 "${BOWTIE2_OUTPUT}/Human_host_removed.read1.fastq.gz" \
    -2 "${BOWTIE2_OUTPUT}/Human_host_removed.read2.fastq.gz" \
    -0 /dev/null -s /dev/null -n

    #echo "SAMPLE_NAME,Host_Contaminant,Total_mapped_reads,Mapped_rate(%),Total_mapped_paired_reads,Mapped_Paired_rate(%)" >> ${LOG_DIR}/Host_Contaminants_stats.csv

    LINE="${SAMPLE_NAME},Human,"
    mapped_reads=`grep "mapped ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d " " -f 1`
    mapped_rate=`grep "mapped ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
    mapped_paired_reads=`grep "paired ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d " " -f 1`
    mapped_paired_rate=`grep "paired ("  "${LOG_DIR}/Human_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
    LINE+="${mapped_reads},${mapped_rate},${mapped_paired_reads},${mapped_paired_rate}"
    echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv
else
    LINE="${SAMPLE_NAME},Human,0,0,0,0"
    echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv
fi

if [ ${alignmentRate} -lt 95 ] || [ ${mouse} -eq 1 ];
then
  # download genome index
  #aws s3 sync --quiet s3://maf-versioned/ninjamap/bowtie2/mouse/ ${DOWNLOAD_DB}/GRCm39/
  # link genome index
  ln -s /mnt/efs/databases/Bowtie2/mouse ${DOWNLOAD_DB}/GRCm39
  #aws s3 cp --quiet s3://czbiohub-microbiome/ReferenceDBs/bowtie2/GRCm39.zip ${DOWNLOAD_DB}/GRCm39.zip
  #unzip "${DOWNLOAD_DB}/GRCm39.zip"
  # bowtie2 mapping against host sequence database, keep both aligned and unaligned reads (paired-end reads)
  bowtie2 -p ${mycpu} -x "${DOWNLOAD_DB}/GRCm39/mouse" -1 "${QC_FASTQ}/read1_trimmed.fastq.gz" -2 "${QC_FASTQ}/read2_trimmed.fastq.gz" -S "${QC_FASTQ}/Mouse_mapped_and_unmapped.sam"
  #convert file .sam to .bam
  samtools view -bS "${QC_FASTQ}/Mouse_mapped_and_unmapped.sam" > "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam"
  samtools flagstat "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam" > "${LOG_DIR}/Mouse_contamination-stat.txt"
  # SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
  samtools view -b -f 12 -F 256 "${QC_FASTQ}/Mouse_mapped_and_unmapped.bam" > "${QC_FASTQ}/Mouse_bothReadsUnmapped.bam"
  # sort bam file by read name (-n) to have paired reads next to each other
  samtools sort -n -@ ${mycpu} "${QC_FASTQ}/Mouse_bothReadsUnmapped.bam" -o "${QC_FASTQ}/Mouse_bothReadsUnmapped_sorted.bam"
  samtools fastq -@ ${mycpu} "${QC_FASTQ}/Mouse_bothReadsUnmapped_sorted.bam" \
  -1 "${BOWTIE2_OUTPUT}/Mouse_host_removed.read1.fastq.gz" \
  -2 "${BOWTIE2_OUTPUT}/Mouse_host_removed.read2.fastq.gz" \
  -0 /dev/null -s /dev/null -n

  LINE="${SAMPLE_NAME},Mouse,"
  mapped_reads=`grep "mapped ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d " " -f 1`
  mapped_rate=`grep "mapped ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
  mapped_paired_reads=`grep "paired ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d " " -f 1`
  mapped_paired_rate=`grep "paired ("  "${LOG_DIR}/Mouse_contamination-stat.txt" | cut -d "%" -f 1 | cut -d"(" -f 2`
  LINE+="${mapped_reads},${mapped_rate},${mapped_paired_reads},${mapped_paired_rate}"
  echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv
else
  LINE="${SAMPLE_NAME},Mouse,0,0,0,0"
  echo "$LINE" >> ${LOG_DIR}/Host_Contaminants_stats.csv
fi


echo "NinjaMap completed."
ls ${LOCAL}
du -sh ${LOCAL}
date
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > ${LOCAL_OUTPUT}/job.complete
echo "Live long and prosper" >> ${LOCAL_OUTPUT}/job.complete
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
