#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154
############################
## NOT WORKING #############
############################


set -e
set -u
set -o pipefail

SAMPLE_NAME=$1
abundance=$2
summary_coverage=$3

cut -f2,3 "${summary_coverage}" | paste "${abundance}" -  | awk '{print $1","$2","$3}' > "${SAMPLE_NAME}.ninjaMap.abundance.csv"

#cut -f2,3 "${summary_coverage}"  > extracted_columns.tsv
#cat extracted_columns.tsv | tr "\t" "," > extracted_columns.csv
#paste -d ',' "${abundance}" extracted_columns.csv > merged_output.csv
#mv merged_output.csv "${SAMPLE_NAME}.ninjaMap.abundance.csv"