#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

# --- Coverage Processing Script ---
# This script handles parallel coverage calculation for singular and escrow BAM files.
# Split from coverage_debug_parallel_2.sh for modularity.

# --- Configuration & Input ---
# These variables must be set in the environment or before calling the script.

# : "${VAR:?Error message}" ensures the script fails if a variable is not set.
: "${in_abundance:?in_abundance is not set}"
: "${singular_bam:?singular_bam is not set}"
: "${escrow_bam:?escrow_bam is not set}"
: "${REF:?REF variable is not set}"
: "${SAMPLE_NAME:?SAMPLE_NAME variable is not set}"

NINJA_OUTPUT="."
GENOME_COV_OUTPUT="."

################################################################################
# PARALLELIZED FUNCTION: Calculate coverage for a single genome
# This function processes one genome and can be run in parallel
# Arguments:
#   $1 - genome name
#   $2 - bed file path
#   $3 - genomes file path
#   $4 - output type (singular/escrow)
################################################################################
calculate_coverage_for_genome() {
    local genome_name=$1
    local bed_file=$2
    local genomes_file=$3
    local output_type=$4
    
    # Use process ID and genome name for unique temporary files
    local tmp_bed="tmp.${genome_name}.$$.bed"
    local tmp_genome="tmp.${genome_name}.$$.genome"
    local tmp_bedg="tmp.${genome_name}.$$.bedg"
    
    # Clean up temporary files on exit
    trap 'rm -f "${tmp_bed}" "${tmp_genome}" "${tmp_bedg}"' RETURN

    # Filter bed file for current genome
    if [ $(grep -c "^$genome_name" "${bed_file}") -eq 0 ]; then
        > "${tmp_bed}"
    else
        grep "^$genome_name" "${bed_file}" | sort | uniq | cut -f1,2,3 > "${tmp_bed}"
    fi
    
    # Get genome info
    grep "^$genome_name" "${genomes_file}" > "${tmp_genome}"
    
    # Run bedtools genomecov
    bedtools genomecov -i "${tmp_bed}" -g "${tmp_genome}" -bga > "${tmp_bedg}"
    
    # Calculate coverage statistics
    zero=$(awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero+0}' "${tmp_bedg}")
    nonzero=$(awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero+0}' "${tmp_bedg}")
    
    if [ -z "$zero" ]; then zero=0; fi
    if [ -z "$nonzero" ]; then nonzero=0; fi
    
    if [ $((zero + nonzero)) -eq 0 ]; then
        cov_pct=0
    else
        cov_pct=$(bc <<< "scale=4; ($nonzero / ($zero + $nonzero))*100")
    fi
    
    totalbases=$(awk '{t+=$4*($3-$2)} END {print t+0}' "${tmp_bedg}")
    if [ $nonzero -eq 0 ] || [ $totalbases -eq 0 ]; then
        cov_depth=0
    else
        cov_depth=$(bc <<< "scale=4; ($totalbases / $nonzero)")
    fi
    
    # Output result
    printf "${genome_name}\t${cov_pct}\t${cov_depth}\t${nonzero}\n"
}

# Export function for use by parallel
export -f calculate_coverage_for_genome

################################################################################
# MAIN COVERAGE PROCESSING
################################################################################

echo "=== Starting Coverage Processing ==="
coreN=$(grep -c ^processor /proc/cpuinfo)
echo "Using ${coreN} cores for parallel processing"

# Prepare reference files
samtools faidx ${REF} && cut -f1,2 ${REF}.fai > Ninja.genomes
cut -f1 Ninja.genomes | sed 's/_Node.*//' | sort | uniq > DBGenomeNameList.txt

echo "Processing $(wc -l < DBGenomeNameList.txt) genomes"

# --- SINGULAR COVERAGE PROCESSING ---
echo "=== Compute Singular coverage and depth ==="
s="singular"
samtools sort -@ ${coreN} ${singular_bam} > ${SAMPLE_NAME}.singular_sorted.bam
bedtools bamtobed -i ${SAMPLE_NAME}.singular_sorted.bam | cut -f1,2,3 > ${SAMPLE_NAME}.${s}.bed
samtools index ${SAMPLE_NAME}.singular_sorted.bam
printf "Strain_Name\tSingular_Coverage\tSingular_Depth\tSingular_Bases\n" > ${s}_summary_depth.tsv

echo "Processing singular genomes in parallel using ${coreN} cores..."

# Check if GNU parallel is available, fallback to xargs if not
if command -v parallel >/dev/null 2>&1; then
    # Use GNU parallel
    cat DBGenomeNameList.txt | parallel -j ${coreN} --bar calculate_coverage_for_genome {} ${SAMPLE_NAME}.${s}.bed Ninja.genomes ${s} >> ${s}_summary_depth.tsv
else
    echo "GNU parallel not found, using xargs with background processes..."
    # Fallback to xargs with parallel processing
    cat DBGenomeNameList.txt | xargs -n 1 -P ${coreN} -I {} bash -c 'calculate_coverage_for_genome "$@"' _ {} ${SAMPLE_NAME}.${s}.bed Ninja.genomes ${s} >> ${s}_summary_depth.tsv
fi

# Merge singular results
merge_abundance_coverage_table.py ${SAMPLE_NAME} ${in_abundance} ${s}_summary_depth.tsv

# --- ESCROW COVERAGE PROCESSING ---
echo "=== Compute Escrow coverage and depth ==="
s="escrow"
samtools sort -@ ${coreN} ${escrow_bam} > ${SAMPLE_NAME}.escrow_sorted.bam
bedtools bamtobed -i ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam | cut -f1,2,3 > ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed
samtools index ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam
printf "Strain_Name\tEscrow_Coverage\tEscrow_Depth\tEscorw_Bases\n" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv

echo "Processing escrow genomes in parallel using ${coreN} cores..."

# Check if GNU parallel is available, fallback to xargs if not
if command -v parallel >/dev/null 2>&1; then
    # Use GNU parallel
    cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt | parallel -j ${coreN} --bar calculate_coverage_for_genome {} ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed ${GENOME_COV_OUTPUT}/Ninja.genomes ${s} >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
else
    echo "GNU parallel not found, using xargs with background processes..."
    # Fallback to xargs with parallel processing
    cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt | xargs -n 1 -P ${coreN} -I {} bash -c 'calculate_coverage_for_genome "$@"' _ {} ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed ${GENOME_COV_OUTPUT}/Ninja.genomes ${s} >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
fi

# Merge escrow results with previous abundance data
merge_abundance_coverage_table.py ${SAMPLE_NAME} ${SAMPLE_NAME}.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${s}_summary_depth.tsv

# Clean up temporary files
rm -f ${GENOME_COV_OUTPUT}/tmp.bed ${GENOME_COV_OUTPUT}/tmp.genome ${GENOME_COV_OUTPUT}/tmps.bedg 2>/dev/null || true

echo "=== Coverage Processing Complete ==="
echo "Results saved in: ${SAMPLE_NAME}.ninjaMap.abundance.csv"
