#!/usr/bin/env nextflow

// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run NinjaMap pipeline against a specific database

    Required Arguments:
      --seedfile       file      a file contains sample name, reads1 and reads2
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
  .fromPath(params.seedfile)
  .ifEmpty { exit 1, "Cannot find the input seedfile" }


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
 	.fromPath(params.seedfile)
 	.ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
  .splitCsv(header: ['sample', 'reads1', 'reads2'], sep: ',', skip: 1)
 	.map{ row -> tuple(row.sample, row.reads1, row.reads2)}
 	.set { seedfile_ch }

seedfile_ch.into { seedfile_ch1; seedfile_ch2 }

/*
 * Run NinjaMap preprocessing
 */
process ninjaMap_preprocessing {

    container params.container
    cpus 16
    memory 64.GB
    publishDir "${output_path}", mode:'copy'

    input:
	  tuple val(sample), val(reads1), val(reads2) from seedfile_ch1

    output:
    file "${sample}" into out_ch

    script:
    """
    export sampleRate="${params.sampleRate}"
    export coreNum="${params.coreNum}"
    export memPerCore="${params.memPerCore}"
    export fastq1="${reads1}"
    export fastq2="${reads2}"
    export REFDBNAME="${params.db_prefix}"
    export S3DBPATH="s3://maf-versioned/ninjamap/Index/${params.db}/db/"
    export S3OUTPUTPATH="${output_path}/${sample}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    /work/ninjaMap_nf_preprocessing.sh
    echo "Done" > ${sample}
    """
}


/*
 * Run NinjaMap
 */
process ninjaMap {

    container params.container
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
    export S3OUTPUTPATH="${output_path}/${f}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    /work/ninjaMap_nf_core.sh
    """
}
