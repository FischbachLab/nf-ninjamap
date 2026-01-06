#!/bin/bash

# version = '1.2.0'
# Generate a seedfile for all samples within an s3 project bucket for ninjaMap 
# File size less than 100MB will be excluded in the seedfile
# Example:
# bash  make_s3_seedfile.sh s3://maf-sequencing/Illumina/240214_A01679_0144_BHWK3TDSX7/MITI-001_DrugSubstance/ MITI-001_DrugSubstance_seedfile.csv 


set -euoE pipefail

S3PATH=${1:?"Specity an s3 path,e.g., s3://genomics-workflow-core/Results/Basespace/NextSeq/20251231_MITI-001-DP-USP-Breakthrough-WGS_H7NK7AFXC/"} 
SEEDFILE=${2:?"Specity a seedfile name, e.g., myproject_seedfile.csv"}
profile=${3:-"default"}

ignored_sample="ctrl"

tmp_sample=$(mktemp /tmp/tmp_sample_list.XXXXXX)

bucket=$(echo "$S3PATH" | cut -d/ -f3)
key=$(echo "$S3PATH" | cut -d/ -f4-)

#echo "$bucket"
#echo "$key"

thres="104857600"

echo -e "sampleName,R1,R2" > $SEEDFILE

echo "Using $profile s3 account"
echo "Excluded samples in $S3PATH"

aws --profile $profile s3api list-objects-v2 --bucket "$bucket" --prefix "$key" --query "Contents[].Key" --output text | tr '\t' '\n' | xargs -n1 basename > $tmp_sample
# Loop through unique sample IDs
for f in $(cat $tmp_sample | grep -E "R1.*fastq.gz|R1.*fq.gz"); do
   size=$(aws --profile $profile s3api head-object --bucket "$bucket" --key "$key${f}" --query 'ContentLength' --output text)
   # check sample size and name
   if [[ "$size" -lt $thres && ! "$f" =~ negctrl ]]; then
        #sampleID=$(echo "$f" | sed -E 's/_S[0-9]+*//')
	sampleID=$(echo "$f" | sed -E 's/_S[0-9]+_R[12](_001)?(\.qcd)?\.(fastq|fq)\.gz//')

        echo "$sampleID,$S3PATH${f} $(awk "BEGIN {printf \"%.2f MB\", $size/1048576}")"
   else

    # Skip if no files match
    #[ -e "$f" ] || continue

    # Extract sample base (remove R1/R2 suffix + extensions)
    sample=$(echo "$f" | sed -E 's/_R[12](_001)?(\.qcd)?\.(fastq|fq)\.gz//')
    #echo $sample
    # Extract suffix pattern (everything after R1/R2)
    suffix=$(echo "$f" | sed -E "s/^${sample}_R[12]//")

    R1="${S3PATH}${sample}_R1${suffix}"
    R2="${S3PATH}${sample}_R2${suffix}"

    # Clean sampleID (remove AJ_date prefix if desired)
    sampleID=$(echo "$sample" | sed -E 's/_S[0-9]+//')

    echo "$sampleID,$R1,$R2" >> $SEEDFILE
   fi
done


# remove tmp file
rm $tmp_sample

# Validate the s3 path
s3_path=`tail -n 1 $SEEDFILE | cut -d, -f2`

if aws s3 --profile $profile ls "$s3_path" > /dev/null 2>&1; then
  # Path is valid, do nothing
  :
else
  # Path is invalid, print and exit
  echo "Invalid S3 path: $s3_path"
  echo "Please check the seedfile!"
  exit 1
fi

