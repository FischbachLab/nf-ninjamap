# NinjaMap - beta

Comes in two flavors.

## Narrow

- Calculate strain abundance, for a given database
- **Use Case** : Manufacturing QC purposes. When you want precise control over who and how much of an expected strain is present in your concoction.

## Broad (roadmap)

- Calculate strain abundance, for a given database and related strains
- **Use Case** : Exploratory purposes. When you wish to find out who and how much of an unknown strain is present in your concoction

## Usage

1. Create and index your database with `ninjaIndex.py`
    1. This tool will accept a directory of your reference genomes (one genome per file) and return a `binmap` file along with a concatenated fasta file of your references. Here are the steps involved in indexing the database.
        1. Run biogrinder on individual genomes to obtain fastq files with predetermined uniform coverage (usually 10x)
            See: the `sunitjain/biogrinder` docker image.
        2. Align the fastq files to the concatenated database individually.
            See: `sunitjain/bowtie2` docker image.
        3. Merge the bam files obtained by individually aligning the fastq files.
            - Download: `aws s3 cp s3://czbiohub-microbiome/Sunit_Jain/Synthetic_Community/ninjaMap/20190730_GroundTruth_uniform100x_NM_MR/ bam_files/ --recursive --exclude '*' --include '*/bowtie2/*bam*'`
            - List: `ls bam_files/*/bowtie2/*.bam > bamfiles.list`
            - Merge: `bamtools merge -list bamfiles.list -out uniform100x.merged.bam &> bamtools.log &`
            - **NOTE:** you might need to install bamtools. Use: `conda install -y -c bioconda bamtools`
        4. Calculate the uniqueness of the genome in the database along with other contigs related metadata.
        5. Concatenate the reference genomes into a single fasta file.
    2. Use this concatenated fasta for all your alignment needs with this database.
    3. This step only needs to be executed once per database.
    4. **TODO**: Convert this to a Nextflow pipeline

## Requirements

1. AWS account with access to:
    - S3
    - AWS Batch with ability to use instances with at least 8 vcpus and 32Gb memory.
2. Local Nextflow setup. (Future release)

## To Do

1. Port as much of pipeline to Nextflow as possible
    - Start with `ninjaIndex`, the current setup has a circular dependency on reference fasta files!
2. Develop a CloudFormation template with appropriate bucket, job queue and compute environment access.
3. Make the entire process more coherent and streamlined.

## Questions / Concerns

- [Sunit Jain](microbiome.ninja) (dev)
- Xiandong Meng (dev)
- Brian Yu (dev)
- Michael Fischbach (PI)
