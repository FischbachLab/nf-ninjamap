#!/bin/bash

#check if a ninjaMap job done by each sample name
# Usage: check_incomplete_ninjaMap_job.sh s3://genomics-workflow-core/Results/Ninjamap/MITI-001-DPv4/MITI-001-DPv4_157dropouts/ ../MITI-001-DPv4-strain_dropout.seedfile.csv

set -euoE pipefail


S3PATH=${1:?"Provide an output s3 path"}
#SAMPLENAME=${2:?"Supply an input sample name file"}
SEEDPATH=${2:?"Supply an input seedfile path"}

while IFS=',' read -ra line; do
  SAMPLE_NAME=${line[0]}
  if ! aws s3 ls ${S3PATH}${SAMPLE_NAME}/job.complete 2>>/dev/null | grep -q "job.complete"; then
 # echo "${line} does not exist or access denied"
 # echo ${line} >> misssed_sample.list
    grep  "${SAMPLE_NAME}" ${SEEDPATH} 
  fi
done<${SEEDPATH}

