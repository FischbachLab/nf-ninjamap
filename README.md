Sample scripts
====================

The simple scripts shows how to run ninjamap jobs under the Nextflow framework.

## Local job example
```{bash}
nextflow run -resume main.nf --reads1 's3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R1_001.fastq.gz' --reads2 's3://czbiohub-microbiome/Xiandong_Meng/test/data/blank_S41_R2_001.fastq.gz' --db  --prefix dtest1 -profile docker
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
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/Ninjamap/biohub" "
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
"--output_path", "s3://genomics-workflow-core/Pipeline_Results/Ninjamap/" "
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
```{bash}
aws batch submit-job \
    --job-name nf-ninjamap \
    --job-queue priority-maf-pipelines \
    --job-definition nextflow-production \
    --container-overrides command="fischbachlab/nf-ninjamap, \
"--seedfile", "s3://genomics-workflow-core/Results/Ninjamap/20221018/seedfile2.csv", \
"--db","HCom2_20221117", \
"--db_prefix", "HCom2", \
"--output_path", "s3://genomics-workflow-core/Results/Ninjamap/HCom2/20221018", \
"--sampleRate", "0.5" "
``

Output files for each sample
====================

The output files are organized into 4 folders.

## bowtie2 folder

The alignment file of all input reads aligned the defined community database in the bam format

## Logs folder

The running logs of various scripts

## ninjaMap folder

1. \*.ninjaMap.abundance.csv: this file shows the statistics of the abundance, coverage and depth of each strain in the defined community
Strain_Name: strain name
Read_Fraction: abundance in the defined community in percentage
Percent_Coverage: average coverage per strain in percentage
Coverage_Depth: average coverage depth

2. \*.ninjaMap.read_stats.csv: this file shows the statistics of input reads
File_Name: sample name
Reads_Aligned: the number of aligned reads
Reads_wPerfect_Aln: the number of aligned reads perfectly
Reads_wSingular_Votes: the number of reads voted as singular
Reads_wEscrowed_Votes: the number of reads voted as escrow
Discarded_Reads_w_Perfect_Aln: the number of discarded perfect reads

3. \*.ninjaMap.strain_stats.csv: this file shows the various statistics of each strains

4. \*.ninjaMap.votes.csv.gz: the statistics of reads voting (singular or escrow)

## Stats folder
adapter_trimming_stats_per_ref.txt: this file shows the statistics of adapter trimming
read_accounting.csv: this file shows the statistics shows the total number of reads, the number of reads after trimming and the number of aligned reads


Aggregated output files for each study
====================

The aggregated output files are organized into 6 files.

1. \*.covDepth.csv: this file shows the average coverage depth per strain by samples
2. \*.host_contaminants.csv: this file shows the detected host contaminants (Human or Mouse) by samples if the unaligned reads is over 5%
3. \*.long.csv: this is the long format of three files (\*.readFraction.csv, \*.covDepth.csv and\*.percCoverage.csv)
4. \*.percCoverage.csv: this file shows the average coverage per strain in percentage by samples
5. \*.reads_stats.csv: this file shows the reads the statistics in read number by samples
6. \*.readFraction.csv: this file shows the abundance in the defined community in percentage by samples
