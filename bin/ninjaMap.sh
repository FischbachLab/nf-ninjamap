#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail


BOWTIE2_BAM=$1
BOWTIE2_BAI=$2

referenceNameFile=${STRAIN_MAP_FILENAME}

STATS_DIR="Stats"
LOG_DIR="Logs" 
NINJA_OUTPUT="ninjaMap"

mkdir -p ${STATS_DIR} ${NINJA_OUTPUT} ${LOG_DIR}

# Exclude Reads Overlapping BED Regions

if [ -z "${mask_bed} " ]; then
  echo "Enable genome masking option"
  bed_option="-bed ${mask_bed}"
else
  bed_option=""
fi

  # check if adding -coverage option
if [ ${coverage} -eq 1 ]; then
  cov_option="-coverage"
else
  cov_option=""
fi

# get bam file size
bamSizeCheck=$(stat -c %s "${BOWTIE2_BAM}")

all_mapped_reads=2
# 100k?
if [ $bamSizeCheck -lt 100000 ]; then
    samtools flagstat "${BOWTIE2_BAM}" > "${LOG_DIR}/Alignment-stat.txt"
    all_mapped_reads=`grep "mapped ("  "${LOG_DIR}/Alignment-stat.txt" | cut -d " " -f 1`
fi

# make sure the number of aligned reads > 1
if [ $all_mapped_reads -gt 1 ]; then
    ninjaMap_parallel.py \
    -bam ${BOWTIE2_BAM} \
    -bin ${referenceNameFile} \
    -outdir ${NINJA_OUTPUT} \
    ${cov_option} \
    ${bed_option} \
    -msv ${singular_vote} \
    -prefix ${OUTPUT_PREFIX} 
else
    mkdir ninjaMap
    touch "ninjaMap/${OUTPUT_PREFIX}.ninjaMap.abundance.csv"
    exit 0
fi