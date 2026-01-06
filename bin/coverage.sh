#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

SAMPLE_NAME=$1
bam=$2
bai=$3
REF=$4

bamParser.py -bam $bam -id 100 -aln_len 100 -out "${SAMPLE_NAME}.filtered.bam"

pileup.sh in=${SAMPLE_NAME}.filtered.bam out=${SAMPLE_NAME}_coverage.txt overwrite=t delcoverage=f 2>${SAMPLE_NAME}_stats.txt

samtools faidx ${REF} && cut -f1,2 ${REF}.fai > Ninja.genomes
cut -f1  Ninja.genomes | sed  's/_Node.*//' | sort | uniq  > DBGenomeNameList.txt
echo -e "Stain_Name\tPercent_Coverage\tCoverage_Depth" > ${SAMPLE_NAME}_summary_coverage.tsv
for i in $(cat DBGenomeNameList.txt)
do
    grep $i ${SAMPLE_NAME}_coverage.txt | awk '{c+=$6; s+=$2*$3;l+=$3}END{print $1"\t"c/l*100"\t"s/l}' >> ${SAMPLE_NAME}_summary_coverage.tsv
done