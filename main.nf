#!/usr/bin/env nextflow
import groovy.json.JsonOutput

nextflow.enable.dsl=2

include { printParams; get_software_versions } from './modules/house_keeping'
include { read_qc; bowtie2_alignment } from './modules/pre_processing'
include { ninjaMap_abundance; ninjaMap_coverage_general; ninjaMap_merge_abundance_table } from './modules/ninjaMap'
include { ninjaMap_contam; ninjaMap_coverage_debug; ninjaMap_coverage; ninjaMap_debug; aggregate_ninjaMap; ninjaMap_complete } from './modules/post_processing'

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
	exit 1, "Missing the database prefix"
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


seedfile_ch  = Channel
    .fromPath(params.seedfile)
    .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
    .splitCsv(header: ['sampleName', 'R1', 'R2'], sep: ',', skip: 1)
    .map{ row -> tuple(row.sampleName, row.R1, row.R2)}
   // .filter { sampleName, R1, R2 -> !sampleName.toLowerCase().contains("negctrl") }

db_ch = Channel
      .fromPath("${params.db_path}/${params.db}/db/bowtie2_index/*")
      .ifEmpty { exit 1, "Cannot find the ninjaMap DB index files: ${params.db}" }
 
fna_ch = Channel
      .fromPath("${params.db_path}/${params.db}/db/*.fna")
      .ifEmpty { exit 1, "Cannot find the ninjaMap DB geneome file: ${params.db}" }

binmap_ch = Channel
      .fromPath("${params.db_path}/${params.db}/db/*.ninjaIndex.binmap.csv")
      .ifEmpty { exit 1, "Cannot find the ninjaMap DB binmap file: ${params.db}" }
  

workflow {

  printParams()
  get_software_versions()

  seedfile_ch | read_qc 
  
  bowtie2_alignment(db_ch.toSortedList(), read_qc.out.qc_fastq_ch.join(seedfile_ch, by: 0))

  // debugging: print the channel
  //bowtie2_alignment.out.bam_ch.view { sample, bam, bai ->
  //  "Sample: $sample | bam: $bam  | bai: $bai"
 // }
  ninjaMap_abundance(bowtie2_alignment.out.bam_ch, binmap_ch.collect())
  ninjaMap_coverage_general(bowtie2_alignment.out.bam_ch, fna_ch.collect())
  ninjaMap_merge_abundance_table(ninjaMap_abundance.out.abundance.join(ninjaMap_coverage_general.out.coverage, by: 0))

  ninjaMap_contam(read_qc.out.qc_fastq_ch.join(bowtie2_alignment.out.bam_ch, by: 0))

  if(params.coverage){
    ninjaMap_coverage(ninjaMap_abundance.out.abundance_coverage.join(ninjaMap_merge_abundance_table.out.abundance, by: 0), fna_ch.collect())
    if(params.debug){
      ninjaMap_debug(ninjaMap_abundance.out.abundance_coverage.join(ninjaMap_merge_abundance_table.out.abundance, by: 0), fna_ch.collect())
    }
    //ninjaMap_coverage_debug(ninjaMap_abundance.out.abundance_coverage.join(ninjaMap_merge_abundance_table.out.abundance, by: 0), fna_ch.collect())
    ninjaMap_complete(ninjaMap_coverage.out.abundance)
  }else{
    ninjaMap_complete(ninjaMap_merge_abundance_table.out.abundance)
  }

  aggregate_ninjaMap(ninjaMap_complete.out.txt.toSortedList())
}


 workflow.onComplete {
     println "Pipeline completed!"
     println "Started at  $workflow.start"
     println "Finished at $workflow.complete"
     println "Time elapsed: $workflow.duration"
     println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
 }
