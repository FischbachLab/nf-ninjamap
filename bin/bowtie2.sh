#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

coreNum="${coreNum:-16}"
coreN="${coreN:-15}"
memPerCore="${memPerCore:-2G}"
maxInsert="${maxInsert:-3000}"
maxAlignments="${maxAlignments:-200}"

DB=$1
fastq1=$2
fastq2=$3
raw1=$4

OUTPUTDIR="tmp"
LOG_DIR="Logs"
STATS_DIR="Stats"

BOWTIE2_OUTPUT="bowtie2"
LOCAL_DB_PATH="reference"
TMP_OUTPUTS="tmp"
OUTPUT_PREFIX="${SAMPLE_NAME}.sortedByCoord"
BOWTIE2_DB="${LOCAL_DB_PATH}/bowtie2_index"
#{REFDBNAME}


mkdir -p ${LOG_DIR} ${BOWTIE2_DB} ${BOWTIE2_OUTPUT} ${LOCAL_DB_PATH} ${TMP_OUTPUTS} ${STATS_DIR}

#BOWTIE2_DB=/bowtie2_index/${REFDBNAME}

cp -r ${DB}/* ${LOCAL_DB_PATH}/bowtie2_index/

# Copy genome reference over  ${params.db_path}/${params.db}/db/
# Check the reference location compgen -G "/path/to/dir/*" > /dev/null
if [ -n "$(ls -A ${BOWTIE2_DB})" ]
then
	echo "Bowtie2 index files are available" > ${LOG_DIR}/read_mapping.log.txt
else
  echo "Copying hCom2 index files from the external URL" >> ${LOG_DIR}/read_mapping.log.txt
  DBPATH="https://zenodo.org/record/7872423/files/hCom2_20221117.ninjaIndex.tar.gz"
  wget $DBPATH -P ${LOCAL_DB_PATH}/
  tar -xzvf ${LOCAL_DB_PATH}/hCom2_20221117.ninjaIndex.tar.gz  -C ${LOCAL_DB_PATH}/
  mv ${LOCAL_DB_PATH}/HCom2_20221117/db/* ${LOCAL_DB_PATH}/
fi

# Try to use all available cores to speed up alignments
mycpu=$(bc <<< "(`grep -c ^processor /proc/cpuinfo` - 2)" )

thread=$(( (mycpu+2) / 2 ))
# bowtie2 alignment returning multiple alignments and using longer max insert size limites
# output samtools bam file with only properly aligned paired reads.
# added(05/30/2024): output paired-end reads that fail to align concordantly.
# --no-overlap \
#  --un-conc-gz
bowtie2 \
    --very-sensitive \
    -X ${maxInsert} \
    -k ${maxAlignments} \
    --threads ${mycpu} \
    -x ${LOCAL_DB_PATH}/bowtie2_index/${REFDBNAME} \
    --no-mixed \
    --no-discordant \
    --end-to-end \
    --no-unal \
    --un-conc-gz ${BOWTIE2_OUTPUT}/${SAMPLE_NAME}_unmapped_include_overlap_R%.fastq.gz \
    -1 ${fastq1} \
    -2 ${fastq2} | \
    samtools view \
        -@ ${thread} \
        -bh \
        -o ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam - |\
    tee -a ${LOG_DIR}/read_mapping.log.txt

samtools sort \
    -n \
    -@ ${thread} \
    -T ${OUTPUTDIR} \
    -O BAM \
    ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam |\
samtools fixmate \
    -O BAM \
    -cm \
    - - | \
samtools sort \
    -@ ${thread} \
    -T ${OUTPUTDIR} \
    -o ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam

# samtools sort -@ ${coreN} -o ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam ${TMP_OUTPUTS}/${SAMPLE_NAME}.bam
samtools index -@ ${mycpu} ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam

# Tabulate read count
#totalReads=$(( $( zcat ${fastq1}| wc -l ) / 4 ))
readsAfterTrim=$(( $( zcat ${fastq2}| wc -l ) / 4 ))

# Not working on duplicated reads using flagstat instead
#uniqueReads=$( samtools view -f 0x40 ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam | cut -f1 | sort -u | wc -l )
uniqueReads=$( samtools flagstat ${BOWTIE2_OUTPUT}/${OUTPUT_PREFIX}.bam | awk 'NR == 7 {print $1}' )

#totalReads="${totalReads%%*( )}"
totalReads=$(( $(zcat "${raw1}" | wc -l ) / 4 ))

echo 'Sample_Name,Total_Fragments,Fragments_After_Trim,Fragments_Aligned' > ${STATS_DIR}/read_accounting.csv
echo ${SAMPLE_NAME}','${totalReads}','${readsAfterTrim}','${uniqueReads} >> ${STATS_DIR}/read_accounting.csv

#just echo the variable, Nextflow captures it from stdout and emits it as a val.
#echo "${uniqueReads}"

#echo ${uniqueReads} > temp.txt