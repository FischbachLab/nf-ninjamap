#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail


RAW_FASTQ="raw_fastq"
QC_FASTQ="qc_fastq"
LOG_DIR="Logs"
STATS_DIR="Stats"

mkdir ${RAW_FASTQ} ${QC_FASTQ} ${LOG_DIR} ${STATS_DIR}

# Constant definitions for bbduk
trimQuality="${trimQuality:-30}"  # was 25
minLength=${minLength:-50}
kmer_value=${kmer_value:-23}
min_kmer_value=${min_kmer_value:-11}
adapterFile="adapters,phix"

echo "Starting to Process Sample: "${SAMPLE_NAME} > ${LOG_DIR}/bbduk.log.txt

totalReads=$(( $(zcat "${fastq1}" | wc -l ) / 4 ))

mv ${fastq1} ${RAW_FASTQ}/read1.fastq.gz
mv ${fastq2} ${RAW_FASTQ}/read2.fastq.gz

##################################################################################################
# check if there are duplicated headers in raw data
##################################################################################################
# randomly sample 5% reads
reformat.sh samplerate=0.05 sampleseed=123 fixheaders=t in=${RAW_FASTQ}/read1.fastq.gz  in2=${RAW_FASTQ}/read2.fastq.gz  out=${RAW_FASTQ}/sampled_10k_R1.fastq.gz out2=${RAW_FASTQ}/sampled_10k_R2.fastq.gz    
   
header1_count=$(zcat ${RAW_FASTQ}/sampled_10k_R1.fastq.gz| awk 'NR%4==1 {print $0}' | wc -l)
header2_count=$(zcat ${RAW_FASTQ}/sampled_10k_R1.fastq.gz | awk 'NR%4==1 {print $0}' | sort | uniq | wc -l)

if [ "${header1_count}" = "${header2_count}" ]; then
  echo "No duplicated read headers found." >> ${LOG_DIR}/bbduk.log.txt
  mv ${RAW_FASTQ}/read1.fastq.gz ${RAW_FASTQ}/deduped_read1.fastq.gz
  mv ${RAW_FASTQ}/read2.fastq.gz ${RAW_FASTQ}/deduped_read2.fastq.gz
else
  echo "There are duplicated read headers. Deduplicating reads ....."
  # fix: dereplicate reads in raw fastq by headers
  zcat ${RAW_FASTQ}/read1.fastq.gz | awk '{if(NR%4==1) { if(!seen[$0]++) {print $0; getline; print $0; getline; print $0; getline; print $0;} } }' > ${RAW_FASTQ}/deduped_read1.fastq
  zcat ${RAW_FASTQ}/read2.fastq.gz | awk '{if(NR%4==1) { if(!seen[$0]++) {print $0; getline; print $0; getline; print $0; getline; print $0;} } }' > ${RAW_FASTQ}/deduped_read2.fastq
  reformat.sh \
  in=${RAW_FASTQ}/deduped_read1.fastq \
  in2=${RAW_FASTQ}/deduped_read2.fastq \
  out=${RAW_FASTQ}/deduped_read1.fastq.gz \
  out2=${RAW_FASTQ}/deduped_read2.fastq.gz
fi

##################################################################################################
# Sampling reads based on the option
##################################################################################################
# Downsample reads to get results faster

if [ "$(echo "$sampleRate == 1" | bc)" -ne 1 ]; then
    echo "***Start downsampling reads($sampleRate)***" 
    reformat.sh \
    samplerate=${sampleRate} \
    sampleseed=1772 \
    in="${RAW_FASTQ}/deduped_read1.fastq.gz" in2="${RAW_FASTQ}/deduped_read2.fastq.gz" \
    out="${RAW_FASTQ}/read1_sampled.fastq.gz" out2="${RAW_FASTQ}/read2_sampled.fastq.gz" 
elif [ "$(echo "$sampleRate == 1" | bc)" -eq 1 ] && [ "$autoSampling" = true ]; then
    #R1_count=$(zcat ${RAW_FASTQ}/deduped_read1.fastq.gz| awk 'NR%4==1 {print $0}' | wc -l)
    R1_size=$(du -shD --block-size=1M "${RAW_FASTQ}/deduped_read1.fastq.gz" | cut -f 1)
    if [ $R1_size -gt 5102 ]; then # 5GB size threshold
        echo "***Start auto-downsampling reads to 50 million reads***" 
        reformat.sh \
        samplereadstarget=50000000 \
        sampleseed=1772 \
        in="${RAW_FASTQ}/deduped_read1.fastq.gz" in2="${RAW_FASTQ}/deduped_read2.fastq.gz" \
        out="${RAW_FASTQ}/read1_sampled.fastq.gz" out2="${RAW_FASTQ}/read2_sampled.fastq.gz" 
    else
        echo "***Use all input reads***"
        mv "${RAW_FASTQ}/deduped_read1.fastq.gz" "${RAW_FASTQ}/read1_sampled.fastq.gz"
        mv "${RAW_FASTQ}/deduped_read2.fastq.gz" "${RAW_FASTQ}/read2_sampled.fastq.gz" 
    fi
else 
     echo "***Use all input reads2***"
    mv "${RAW_FASTQ}/deduped_read1.fastq.gz" "${RAW_FASTQ}/read1_sampled.fastq.gz"
    mv "${RAW_FASTQ}/deduped_read2.fastq.gz" "${RAW_FASTQ}/read2_sampled.fastq.gz" 
fi


# Use bbduk to trim reads, -eoom exits when out of memory

bbduk.sh -Xmx12g -eoom \
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



#echo "${totalReads}" 



