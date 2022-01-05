#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run NinjaMap pipeline against a specific database

    Required Arguments:
      --reads1        R1        Forward reads file path ( paired-end library )
      --reads2        R2        Reverse reads file path ( paired-end library )
      --db            db_name   NinjaMap database name
      --db_prefix     db_prefix NinjaMap database prefix
      --output_path   path      Output s3 path

    Options:
      --sampleRate    num   Sampling rate (0-1)
      --coreNum       num   Number of cores (e.g. 15)
      --memPerCore    num   Memory per core (e.g., 1G)
      -profile        docker run locally


    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Set default options
if (params.db_prefix == "null") {
	exit 1, "Missing a database prefix"
}

if (params.db == "null") {
	exit 1, "Missing the a database name"
}

Channel
  .fromPath(params.reads1)
  .ifEmpty { exit 1, "Cannot find fastq R1 file" }

Channel
  .fromPath(params.reads2)
  .ifEmpty { exit 1, "Cannot find fastq R2 file" }



/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

//core = params.
//mem =
//srate=

def output_path = "${params.output_path}"
//def output_path = "s3://genomics-workflow-core/Pipeline_Results/NinjaMap/${params.output_prefix}"

//println output_path
/*
 * Given the query parameter creates a channel emitting the query fasta file(s),
 * the file is split in chunks containing as many sequences as defined by the parameter 'chunksize'.
 * Finally assign the result channel to the variable 'fasta_ch'
 */
Channel
    .fromPath(params.reads1)
    .set { read1_ch }

Channel
    .fromPath(params.reads2)
    .set { read2_ch }

/*
 * Run NinjaMap preprocessing
 */
process ninjaMap_preprocessing {

    container "fischbachlab/nf-ninjamap:latest"
    cpus 16
    memory 64.GB
    publishDir "${output_path}", mode:'copy'

    input:
    file read1 from read1_ch
    file read2 from read2_ch

    output:
    file 'flag.txt' into out_ch
    //file 'tmp_*/Sync/bowtie2/*.bam' into bam_ch
    //file 'tmp_*/Sync/bowtie2/*.bai' into bai_ch

    script:
    """
    export sampleRate="${params.sampleRate}"
    export coreNum="${params.coreNum}"
    export memPerCore="${params.memPerCore}"
    export fastq1="${params.reads1}"
    export fastq2="${params.reads2}"
    export REFDBNAME="${params.db_prefix}"
    export S3DBPATH="s3://maf-versioned/ninjamap/Index/${params.db}/db/"
    export S3OUTPUTPATH="${output_path}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    /work/ninjaMap_nf_preprocessing.sh
    echo "Done" > flag.txt
    """
}


/*
 * Run NinjaMap
 */
process ninjaMap {

    container "fischbachlab/nf-ninjamap:latest"
    cpus 32
    memory { 256.GB * task.attempt }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    //publishDir "${output_path}", mode:'copy'

    input:
    file f from out_ch
    output:

    script:
    """
    export REFDBNAME="${params.db_prefix}"
    export S3DBPATH="s3://maf-versioned/ninjamap/Index/${params.db}/db/"
    export S3OUTPUTPATH="${output_path}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    /work/ninjaMap_nf_core.sh
    """
}
