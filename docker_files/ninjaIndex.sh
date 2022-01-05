#!/bin/bash -x
set -e
set -u
set -o pipefail

START_TIME=$SECONDS
export PATH="/opt/conda/bin:${PATH}"

# INPUTS
# PREFIX="uniform10x"
# S3INPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190720_00_NinjaIndex/${PREFIX}/setup
# S3OUTPUTPATH=s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190720_00_NinjaIndex/${PREFIX}/index

# PREFIX="uniform100x"
# S3INPUTPATH="s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190731_00_NinjaIndex/${PREFIX}/setup"
# S3OUTPUTPATH="s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190731_00_NinjaIndex/${PREFIX}/index"

coreNum="${coreNum:-15}"
echo "${PATH}"
LOCAL=$(pwd)
scriptFolder="./scripts"

# Setup directory structure
OUTPUTDIR=${LOCAL}/tmp_$( date +"%Y%m%d_%H%M%S" )

LOCAL_OUTPUT="${OUTPUTDIR}/Sync"
LOG_DIR="${LOCAL_OUTPUT}/Logs"
NINJA_OUTPUT="${LOCAL_OUTPUT}/ninjaIndex"

S3OUTPUTPATH=${S3OUTPUTPATH%/}
S3INPUTPATH=${S3INPUTPATH%/}
LOCAL_DB_PATH="${OUTPUTDIR}/reference_fasta"

mkdir -p "${OUTPUTDIR}" "${LOCAL_OUTPUT}" "${LOG_DIR}"
mkdir -p "${LOCAL_DB_PATH}" "${NINJA_OUTPUT}"

trap '{aws s3 sync "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}";
    rm -rf ${OUTPUTDIR} ; 
    exit 255; }' 1 

aws s3 sync --quiet "${S3INPUTPATH}" "${OUTPUTDIR}"

python -u "${scriptFolder}/ninjaIndex.py" \
    -bam "${OUTPUTDIR}/${PREFIX}.merged.bam" \
    -fastadir "${LOCAL_DB_PATH}" \
    -prefix "${NINJA_OUTPUT}/${PREFIX}" | tee -a "${LOG_DIR}/${PREFIX}_ninjaIndex.binmap.log"

ls "${LOCAL}"
du -sh "${LOCAL}"
date
######################### HOUSEKEEPING #############################
DURATION=$((SECONDS - START_TIME))
hrs=$(( DURATION/3600 )); mins=$(( (DURATION-hrs*3600)/60)); secs=$(( DURATION-hrs*3600-mins*60 ))
printf 'This AWSome pipeline took: %02d:%02d:%02d\n' $hrs $mins $secs > "${LOCAL_OUTPUT}/job.complete"
echo "Live long and prosper" >> "${LOCAL_OUTPUT}/job.complete"
############################ PEACE! ################################
## Sync output
aws s3 sync --quiet "${LOCAL_OUTPUT}" "${S3OUTPUTPATH}"
