#!/usr/bin/env nextflow
import groovy.json.JsonOutput

nextflow.enable.dsl=2
// If the user uses the --help flag, print the help text below
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run NinjaMap pipeline against a specific database

    Required Arguments:
      --seedfile      file      a file contains sample name, reads1 and reads2
      --db            db_name   NinjaMap database name
      --db_prefix     db_prefix NinjaMap database prefix
      --db_path       db_path   NinjaMap database path
      --output_path   path      Output s3 path

    Options:
      --sampleRate    num   Sampling rate (0-1)
      --coreNum       num   Number of cores (e.g. 15)
      --memPerCore    num   Memory per core (e.g., 2G)
      --coverage      num   Outputting singular & escrow coverage and depth (0 or 1, default 0)
      --minQuality    num   Regions with average quality BELOW this will be trimmed
      --minLength     num   Reads shorter than this after trimming will be discarded
      --debug         num   enable debug mode to generate singular bam files (0 or 1, default 0)
      -profile        docker  run locally



    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

if (params.output_path == "null") {
	exit 1, "Missing the output path"
}

if (params.db == "null") {
	exit 1, "Missing the a database name"
}

if (params.db_prefix == "null") {
	exit 1, "Missing the a database prefix"
}

// coverage must be enabled if enabling the debug mode
if (params.debug && !params.coverage ) {
	  exit 1, "The coverage option must be set to 1 if enabling the debug option"
}

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

def output_path = "${params.output_path}"
//def output_path=s3://genomics-workflow-core/Pipeline_Results/NinjaMap/${params.output_prefix}"

/*
 * Run NinjaMap  16 128 fischbachlab/nf-ninjamap:20220726210531
 */
process ninjaMap {
    tag "$sample"
    container  params.container
    cpus { 32 * task.attempt }
    memory { 250.GB * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${output_path}", mode:'copy'

    input:
    tuple val(sample), val(reads1), val(reads2)

    output:

    script:
    """
    export sampleRate="${params.sampleRate}"
    export coreNum="${params.coreNum}"
    export memPerCore="${params.memPerCore}"
    export fastq1="${reads1}"
    export fastq2="${reads2}"
    export coverage="${params.coverage}"
    export REFDBNAME="${params.db_prefix}"
    export S3DBPATH="${params.db_path}/${params.db}/db/"
    export S3OUTPUTPATH="${output_path}/${sample}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    export trimQuality="${params.minQuality}"
    export minLength="${params.minLength}"
    export debug="${params.debug}"
    export ref_db="${params.ref_db_path}"
    ninjaMap.sh
    """
}

/*
  Save all parameters to a jason file
  https://github.com/nextflow-io/nextflow/discussions/2892
  //echo ${params} > ninjamap-parameters.txt
*/
process printParams {

    container  params.container
    publishDir "${params.output_path}"

    errorStrategy = 'ignore'

    output:
    path "parameters.json"

    script:
    """
    touch parameters.json
    echo '${JsonOutput.toJson(params)}' > parameters.json
    """
}

/*
 * Parse software version numbers
 */
process get_software_versions {

    container  params.container
    errorStrategy 'ignore'
    publishDir "${params.output_path}/pipeline_info"
    //saveAs: {filename ->
    //    if (filename.indexOf(".csv") > 0) filename
    //    else null
    //}

    output:
    path 'software_versions_ninjamap.yaml'
    //path 'software_versions_ninja.yaml'
    //path "software_versions.csv"
    //path "*.txt"

    script:
    // Get all tools to print their version number here
    // scrape_software_versions.py &> software_versions_ninjamap.yaml
     //echo $workflow.manifest.version > v_pipeline.txt
    //echo $workflow.nextflow.version > v_nextflow.txt
    //python --version > v_python.txt
    //printf "pipeline_hash: %s\n" ${workflow.scriptId} >> software_versions_ninjamap.yaml
    """
    printf "nextflow_version: %s\n" ${workflow.nextflow.version} > software_versions_ninjamap.yaml
    printf "pipeline_version: %s\n" ${workflow.manifest.version} >> software_versions_ninjamap.yaml
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}') >> software_versions_ninjamap.yaml
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}') >> software_versions_ninjamap.yaml
    printf "bowtie2_version: %s\n" \$(bowtie2 --version | grep -a bowtie2-align-s | awk '{print \$NF}') >> software_versions_ninjamap.yaml
    printf "python_version: %s\n" \$(python --version | awk '{print \$NF}') >> software_versions_ninjamap.yaml
    """
}




Channel
 .fromPath(params.seedfile)
 .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
 .splitCsv(header: ['sample', 'reads1', 'reads2'], sep: ',', skip: 1)
 .map{ row -> tuple(row.sample, row.reads1, row.reads2)}
 .set { seedfile_ch }

workflow {
  get_software_versions()
  seedfile_ch |  ninjaMap
  printParams()
}
