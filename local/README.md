A sample script shows how to run ninjamap jobs under the Nextflow framework on premise servers.
====================

1. Install [nextflow](https://www.nextflow.io/)
2. Run the following command

## [ninjaMap workflow]

## Command line example for a single sample for files stored on an premise servers using a seedfile.
### Note that you must specify an absolute path as the output_path rather than a relative path.
```{bash}
nextflow run main.nf --seedfile /path/to/test_seedfile.csv --output_path /absolute/path/to/local/output
```

## Test dataset
The test dataset is available in the data directory
```{bash}
data/read1.fastq.gz
data/read2.fastq.gz
```