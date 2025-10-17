#!/usr/bin/env nextflow
import groovy.json.JsonOutput

nextflow.enable.dsl=2

include { ninjaMap; aggregate_ninjaMap;} from './modules/ninjamap'
include { printParams; get_software_versions; } from './modules/house_keeping'

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
      --project       project_name  a project name

    Options:
      --sampleRate    num   Sampling rate (0-1)
      --coreNum       num   Number of cores (e.g. 15)
      --memPerCore    num   Memory per core (e.g., 2G)
      --coverage      num   Outputting singular & escrow coverage and depth (0 or 1, default 0)
      --minQuality    num   Regions with average quality BELOW this will be trimmed
      --minLength     num   Reads shorter than this after trimming will be discarded
      --min_singular_vote num   Minimum singular vote per strain
      --debug         num   Enable debug mode to generate singular bam files (0 or 1, default 0)
      --mask          file  Masked genome regions excluding for alingments in bed format
      -profile        docker  Run locally



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

//def output_path = "${params.output_path}"
//def output_path=s3://genomics-workflow-core/Pipeline_Results/NinjaMap/${params.output_prefix}"


Channel
 .fromPath(params.seedfile)
 .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
 .splitCsv(header: ['sample', 'reads1', 'reads2'], sep: ',', skip: 1)
 .map{ row -> tuple(row.sample, row.reads1, row.reads2)}
 .set { seedfile_ch }

workflow {
  printParams()
  get_software_versions()
  seedfile_ch |  ninjaMap
  aggregate_ninjaMap(ninjaMap.out.toSortedList())
  //aggregate_ninjaMap(ninjaMap.out.abund_ch.concat( ninjaMap.out.stats_ch, ninjaMap.out.contam_ch, ninjaMap.out.read_ch ).toSortedList())
}
