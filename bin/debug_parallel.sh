#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

# --- Debug Processing Script ---
# This script handles parallel debug file processing for singular and escrow data.
# Split from coverage_debug_parallel_2.sh for modularity.

# --- Configuration & Input ---
# These variables must be set in the environment or before calling the script.

: "${singular_bam:?singular_bam is not set}"
: "${escrow_bam:?escrow_bam is not set}"
: "${REF:?REF variable is not set}"
: "${SAMPLE_NAME:?SAMPLE_NAME variable is not set}"


# Optional variables with defaults
bedN="${bedN:-100000}" # bedfile threshold 0.1M 

NINJA_OUTPUT="."
GENOME_COV_OUTPUT="."

################################################################################
# PARALLELIZED FUNCTION: Process debug files for a single genome
# Arguments:
#   $1 - processing type (singular/escrow)
#   $2 - genome name
################################################################################
process_debug_genome() {
    local s=$1
    local genome_name=$2
    
    # Use exported variables instead of local ones to ensure they work in parallel
    local debug_dir="${NINJA_OUTPUT}/debug"
    local coreN=$(grep -c ^processor /proc/cpuinfo)
    
    # Ensure the output directory exists (this needs to be done for each parallel process)
    mkdir -p "${debug_dir}/${s}"
    
    local bed_line=$(grep -c "^$genome_name" "${debug_dir}/all/${SAMPLE_NAME}.${s}.bed" 2>/dev/null || echo "0")
    
    # Only create BED files for genomes that have reads (bed_line > 0)
    if [ ${bed_line} -gt 0 ]; then
        # Create BED file with reads for this genome
        grep "^$genome_name" "${debug_dir}/all/${SAMPLE_NAME}.${s}.bed" | sort | uniq | cut -f1,2,3 > "${debug_dir}/${s}/${genome_name}.bed"
        
        # Only create BAM files if reads are below threshold
        if [ ${bed_line} -lt ${bedN} ]; then
            bedtools intersect -abam "${debug_dir}/all/${SAMPLE_NAME}.${s}_sorted.bam" -b "${debug_dir}/${s}/${genome_name}.bed" > "${debug_dir}/${s}/${genome_name}_${s}.bam"
            samtools index -@ ${coreN} "${debug_dir}/${s}/${genome_name}_${s}.bam"
            echo "Processed ${genome_name} for ${s}: ${bed_line} reads -> BED and BAM files created"
        else
            echo "Processed ${genome_name} for ${s}: ${bed_line} reads -> BED file created (>= ${bedN} threshold, BAM file skipped)"
        fi
    else
        echo "Skipped ${genome_name} for ${s}: 0 reads (no BED file created)"
    fi
}

# Export function and necessary variables for use by parallel
export -f process_debug_genome
export NINJA_OUTPUT
export GENOME_COV_OUTPUT 
export SAMPLE_NAME
export bedN

################################################################################
# MAIN DEBUG PROCESSING
################################################################################

echo "=== Starting Debug Processing ==="
coreN=$(grep -c ^processor /proc/cpuinfo)
echo "Using ${coreN} cores for parallel processing"
echo "BED line threshold: ${bedN}"

samtools faidx ${REF} && cut -f1,2 ${REF}.fai > Ninja.genomes
cut -f1  Ninja.genomes | sed  's/_Node.*//' | sort | uniq  > DBGenomeNameList.txt

echo "Compute Singular..."
s="singular"
samtools sort -@ ${coreN} ${singular_bam} > ${SAMPLE_NAME}.singular_sorted.bam
bedtools bamtobed -i ${SAMPLE_NAME}.singular_sorted.bam | cut -f1,2,3 > ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed
samtools index  ${SAMPLE_NAME}.singular_sorted.bam

echo "Compute Escrow..."
s="escrow"
samtools sort -@ ${coreN} ${escrow_bam} > ${SAMPLE_NAME}.escrow_sorted.bam
bedtools bamtobed -i ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam  | cut -f1,2,3 > ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed
samtools index ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam 

# Check if required files exist
if [ ! -f "${GENOME_COV_OUTPUT}/DBGenomeNameList.txt" ]; then
    echo "ERROR: DBGenomeNameList.txt not found in ${GENOME_COV_OUTPUT}/"
    echo "This file should be created by running the coverage script first."
    exit 1
fi

# Create debug directory structure
echo "Setting up debug directories..."
mkdir -p "${NINJA_OUTPUT}/debug/all"

# Move BAM and BED files to debug/all directory
echo "Moving files to debug directory..."

# Move sorted BAM files (with error handling)
if ls ${NINJA_OUTPUT}/*_sorted.bam 1> /dev/null 2>&1; then
    mv ${NINJA_OUTPUT}/*_sorted.bam "${NINJA_OUTPUT}/debug/all/"
    echo "Moved sorted BAM files"
else
    echo "WARNING: No sorted BAM files found to move"
fi

# Move BAM indices (with error handling)
if ls ${NINJA_OUTPUT}/*_sorted.bam.bai 1> /dev/null 2>&1; then
    mv ${NINJA_OUTPUT}/*_sorted.bam.bai "${NINJA_OUTPUT}/debug/all/"
    echo "Moved BAM index files"
else
    echo "WARNING: No BAM index files found to move"
fi

# Move BED files (with error handling)
if ls ${NINJA_OUTPUT}/*.bed 1> /dev/null 2>&1; then
    mv ${NINJA_OUTPUT}/*.bed "${NINJA_OUTPUT}/debug/all/"
    echo "Moved BED files"
else
    echo "WARNING: No BED files found to move"
fi

# Process debug files for each type (singular, escrow)
for s in {singular,escrow}; do
    echo "=== Processing debug files for ${s} genomes ==="
    
    # Check if the required BED file exists
    if [ ! -f "${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}.bed" ]; then
        echo "WARNING: ${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}.bed not found, skipping ${s} processing"
        continue
    fi
    
    # Check if the required BAM file exists
    if [ ! -f "${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}_sorted.bam" ]; then
        echo "WARNING: ${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}_sorted.bam not found, skipping ${s} processing"
        continue
    fi
    
    mkdir -p "${NINJA_OUTPUT}/debug/${s}"
    
    echo "Processing $(wc -l < ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt) genomes for ${s} in parallel using ${coreN} cores..."
    
    # Check if GNU parallel is available, fallback to xargs if not
    if command -v parallel >/dev/null 2>&1; then
        # Use GNU parallel
        cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt | parallel -j ${coreN} --bar process_debug_genome ${s} {}
    else
        echo "GNU parallel not found, using xargs with background processes..."
        # Fallback to xargs with parallel processing
        cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt | xargs -n 1 -P ${coreN} -I {} bash -c 'process_debug_genome "$@"' _ ${s} {}
    fi
    
    # Report results
    bed_files_created=$(find "${NINJA_OUTPUT}/debug/${s}/" -name "*.bed" | wc -l 2>/dev/null || echo "0")
    bam_files_created=$(find "${NINJA_OUTPUT}/debug/${s}/" -name "*_${s}.bam" | wc -l 2>/dev/null || echo "0")
    echo "Created ${bed_files_created} BED files and ${bam_files_created} BAM files for ${s}"
done

echo "=== Debug Processing Complete ==="
echo "Debug files organized in: ${NINJA_OUTPUT}/debug/"
echo "Directory structure:"
echo "  debug/all/          - Original sorted BAM and BED files"
echo "  debug/singular/     - Per-genome files for singular reads"
echo "  debug/escrow/       - Per-genome files for escrow reads"
