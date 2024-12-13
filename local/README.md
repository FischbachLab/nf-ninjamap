A sample script shows how to run ninjamap jobs under the Nextflow framework on a local computer.
====================

1. Install [nextflow](https://www.nextflow.io/)
2. Run the following command

## [ninjaMap workflow]

## Command line example for a single sample for files stored on a local server using a seedfile
```{bash}
nextflow run main.nf --seedfile /path/to/test_seedfile.csv --output_path /path/to/local/output
```

## Test dataset
The test dataset is available in the data directory
```{bash}
data/read1.fastq.gz
data/read2.fastq.gz
```