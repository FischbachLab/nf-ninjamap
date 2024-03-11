#!/bin/bash

# Generate a seedfile for all samples within an s3 project bucket for ninjaMap 
# Example:
# bash  make_s3_seedfile.sh s3://maf-sequencing/Illumina/240214_A01679_0144_BHWK3TDSX7/MITI-001_DrugSubstance/ MITI-001_DrugSubstance_seedfile.csv 


set -euoE pipefail

S3PATH=${1:?"Specity an s3 path,e.g., s3://maf-sequencing/Illumina/240214_A01679_0144_BHWK3TDSX7/MITI-001_DrugSubstance/"} 
SEEDFILE=${2:?"Specity a seedfile name, e.g., myproject_seedfile.csv"}
profile=${3:-"default"}


sample=$(mktemp /tmp/tmp_sample_list.XXXXXX)

echo -e "sampleName,R1,R2" > $SEEDFILE

if [ $profile = "biohub" ] ; then
   sample1=`aws s3 --profile $profile ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | head -n 1`
    if [[ $sample1 == *"_R1.qcd.fq.gz" ]]; then
        aws s3 --profile $profile ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1.qcd.fq.gz//' > $sample
        for i in $(cat $sample); do
          echo -en "$i," >> $SEEDFILE
          echo "$S3PATH${i}_R1.qcd.fq.gz,$S3PATH${i}_R2.qcd.fq.gz" >> $SEEDFILE
        done
    elif  [[ $sample1 == *"_R1.fastq.gz" ]]; then
        aws s3 --profile $profile ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1.fastq.gz//' > $sample
        for i in $(cat $sample); do
          echo -en "$i," >> $SEEDFILE
          echo "$S3PATH${i}_R1.fastq.gz,$S3PATH${i}_R2.fastq.gz" >> $SEEDFILE
        done
    elif [[ $sample1 == *"_R1_001.fastq.gz" ]]; then
       aws s3 --profile $profile ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1_001.fastq.gz//' > $sample
       for i in $(cat $sample); do
         echo -en "$i," >> $SEEDFILE
         echo "$S3PATH${i}_R1_001.fastq.gz,$S3PATH${i}_R2_001.fastq.gz" >> $SEEDFILE
       done
   fi

else
    sample1=`aws s3 ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | head -n 1` 
    if [[ $sample1 == *"_R1.qcd.fq.gz" ]]; then
        aws s3 ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1.qcd.fq.gz//' > $sample
        for i in $(cat $sample); do
          echo -en "$i," >> $SEEDFILE
          echo "$S3PATH${i}_R1.qcd.fq.gz,$S3PATH${i}_R2.qcd.fq.gz" >> $SEEDFILE
        done
    elif  [[ $sample1 == *"_R1.fastq.gz" ]]; then
        aws s3 ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1.fastq.gz//' > $sample
        for i in $(cat $sample); do
          echo -en "$i," >> $SEEDFILE
          echo "$S3PATH${i}_R1.fastq.gz,$S3PATH${i}_R2.fastq.gz" >> $SEEDFILE  
        done
    elif [[ $sample1 == *"_R1_001.fastq.gz" ]]; then
       aws s3 ls $S3PATH | rev | cut -d' ' -f1 | rev | grep "R1" | sed -e 's/S*_R1_001.fastq.gz//' > $sample
       for i in $(cat $sample); do
         echo -en "$i," >> $SEEDFILE
         echo "$S3PATH${i}_R1_001.fastq.gz,$S3PATH${i}_R2_001.fastq.gz" >> $SEEDFILE
       done
   fi

fi

rm $sample 
