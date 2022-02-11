Hello world script
====================

A simple script showing the ninjamap example for the Nextflow framework.
Output directory is at s3://

```{bash}
nextflow run -resume main.nf --reads1 's3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R1_001.fastq.gz' --reads2 's3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R2_001.fastq.gz' --db  --prefix dtest1 -profile docker
```

```{bash}
aws batch submit-job \
  --job-name nf-ninjamap \
  --job-queue priority-maf-pipelines \
  --job-definition nextflow-production \
  --container-overrides command="s3://nextflow-pipelines/nf-ninjamap, \
"--reads1","s3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R1_001.fastq.gz", \
"--reads2","s3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R2_001.fastq.gz", \
"--db","SCv2_6_20210518", \
"--db_prefix", "SCv2_6", \
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/biohub" "
```


```{bash}
aws batch submit-job \
    --job-name nf-ninjamap \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="s3://nextflow-pipelines/nf-ninjamap, \
"--reads1","s3://nextflow-pipelines/nf-ninjamap/data/read1.fastq.gz", \
"--reads2","s3://nextflow-pipelines/nf-ninjamap/data/read2.fastq.gz", \
"--db","SCv2_6_20210518", \
"--db_prefix", "SCv2_6", \
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/2stages" "
```

# Updated version with a seedfile as the input file
```{bash}
aws batch submit-job \
    --job-name nf-ninjamap \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="s3://nextflow-pipelines/nf-ninjamap, \
"--seedfile", "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/seedfile2.csv", \
"--db","SCv2_6_20210518", \
"--db_prefix", "SCv2_6", \
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/2samples", \
"--sampleRate", "0.1"  "
```
