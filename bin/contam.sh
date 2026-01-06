#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

bam=$1
bai=$2

mycpu=$(grep -c ^processor /proc/cpuinfo | bc)

LOG_DIR="Logs"
DOWNLOAD_DB="DB"

mkdir -p ${LOG_DIR} ${DOWNLOAD_DB}

#echo "Totalreads:"
#echo "${uniqueReads}" 

readsAfterTrim=$(( $( zcat ${read1}| wc -l ) / 4 ))
#echo "${readsAfterTrim}"


#uniqueReads=$(($(samtools view -c -F 0x900 ${bam})/2))
uniqueReads=$( samtools flagstat ${bam} | awk 'NR == 7 {print $1}' )
echo "uReads: ${uniqueReads}"

alignmentRate=$(( uniqueReads * 100 / readsAfterTrim ))
echo "alignmentRate: ${alignmentRate}"

echo "SAMPLE_NAME,Host_Contaminant,Total_mapped_reads,Mapped_rate(%),Total_mapped_paired_reads,Mapped_Paired_rate(%)" > ${LOG_DIR}/Host_Contaminants_stats.csv

    if [ ${alignmentRate} -lt $thres ];
    then
    
        ln -s ${ref}/human ${DOWNLOAD_DB}/GRCh38
        
        bowtie2 -p ${mycpu} -x "${DOWNLOAD_DB}/GRCh38/human" -1 "${read1}" -2 "${read2}" -S "Human_mapped_and_unmapped.sam"
        #bowtie2 -p ${mycpu} -x "${ref}/human/human" -1 "${read1}" -2 "${read2}" -S "Human_mapped_and_unmapped.sam"

        samtools view -bS "Human_mapped_and_unmapped.sam" > "Human_mapped_and_unmapped.bam"
        samtools flagstat "Human_mapped_and_unmapped.bam" > "${LOG_DIR}/Human_contamination-stat.txt"
        
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

    if [ ${alignmentRate} -lt $thres ];
    then
        ln -s ${ref}/mouse ${DOWNLOAD_DB}/GRCm39
        
        bowtie2 -p ${mycpu} -x "${DOWNLOAD_DB}/GRCm39/mouse" -1 "${read1}" -2 "${read2}" -S "Mouse_mapped_and_unmapped.sam"
        #bowtie2 -p ${mycpu} -x "${ref}/mouse/mouse" -1 "${read1}" -2 "${read2}" -S "Mouse_mapped_and_unmapped.sam"

        samtools view -bS "Mouse_mapped_and_unmapped.sam" > "Mouse_mapped_and_unmapped.bam"
        samtools flagstat "Mouse_mapped_and_unmapped.bam" > "${LOG_DIR}/Mouse_contamination-stat.txt"

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