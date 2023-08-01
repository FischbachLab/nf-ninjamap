Sample scripts
====================

The simple scripts shows how to run ninjamap jobs under the Nextflow framework.

## Local job example
```{bash}
nextflow run  main-local.nf --reads1 's3://nextflow-pipelines/nf-ninjamap/data/read1.fastq.gz' --reads2 's3://nextflow-pipelines/nf-ninjamap/data/read2.fastq.gz' --db HCom2_20221117 --db_prefix HCom2 --output_path 's3://Results/ninjamap/local' -work-dir 's3://Results/Ninjamap/local/work_dir' -profile docker
```

## Example 1: aws batch job

```{bash}
aws batch submit-job \
  --job-name nf-ninjamap \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="fischbachlab/nf-ninjamap, \
"--reads1","s3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R1_001.fastq.gz", \
"--reads2","s3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R2_001.fastq.gz", \
"--db","SCv2_4_20210212", \
"--db_prefix", "SCv2_4", \
"--db_path", "s3://maf-versioned/ninjamap/Index", \
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/Ninjamap/db_SCv2_4/biohub" "
```

## Example 2: aws batch job
```{bash}
aws batch submit-job \
    --job-name nf-ninjamap \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="fischbachlab/nf-ninjamap, \
"--reads1","s3://nextflow-pipelines/nf-ninjamap/data/read1.fastq.gz", \
"--reads2","s3://nextflow-pipelines/nf-ninjamap/data/read2.fastq.gz", \
"--db","SCv2_6_20210518", \
"--db_prefix", "SCv2_6", \
"--db_path", "s3://maf-versioned/ninjamap/Index", \
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/Ninjamap/db_SCv2_6/project1" "
```


## Seedfile example
### Note that the seedfile is a CSV (comma-separated values) file with header
### The format of the seedfile is sample_name,short_R1,short_R2

```{bash}
sampleName,R1,R2
Plate1_MITI-001-Mouse_A10_W8_6-1_S394,s3://maf-sequencing/Illumina/221213_A01679_0069_BHLLVHDSX5/Allison_Weakley/MITI-001-BackfillAnalysis/Plate1_MITI-001-Mouse_A10_W8_6-1_S394_R1.fastq.gz,s3://maf-sequencing/Illumina/221213_A01679_0069_BHLLVHDSX5/Allison_Weakley/MITI-001-BackfillAnalysis/Plate1_MITI-001-Mouse_A10_W8_6-1_S394_R2.fastq.gz
```

## Example 3: aws batch job using a seedfile
### Updated version with a seedfile as the input file
### Example full db path: s3://maf-versioned/ninjamap/Index/HCom2_20221117/db/
```{bash}
aws batch submit-job \
    --job-name nf-ninjamap \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="fischbachlab/nf-ninjamap, \
"--seedfile", "s3://genomics-workflow-core/Results/Ninjamap/20221018/seedfile2.csv", \
"--db","HCom2_20221117", \
"--db_prefix", "HCom2", \
"--db_path", "s3://maf-versioned/ninjamap/Index", \
"--output_path", "s3://genomics-workflow-core/Results/Ninjamap/HCom2/20221018", \
"--sampleRate", "0.5" "
```
