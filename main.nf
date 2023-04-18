#!/usr/bin/env nextflow
nextflow.enable.dsl=1
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

/*
Channel
  .fromPath(params.reads1)
  .ifEmpty { exit 1, "Cannot find fastq R1 file" }

Channel
  .fromPath(params.reads2)
  .ifEmpty { exit 1, "Cannot find fastq R2 file" }


 * Defines the pipeline inputs parameters (giving a default value for each for them)
 * Each of the following parameters can be specified as command line options
 */

def output_path = "${params.output_path}"
//def output_path=s3://genomics-workflow-core/Pipeline_Results/NinjaMap/${params.output_prefix}"

//println output_path
/*
Channel
    .fromPath(params.reads1)
    .set { read1_ch }

Channel
    .fromPath(params.reads2)
    .set { read2_ch }
*/

Channel
 .fromPath(params.seedfile)
 .ifEmpty { exit 1, "Cannot find any seed file matching: ${params.seedfile}." }
 .splitCsv(header: ['sample', 'reads1', 'reads2'], sep: ',', skip: 1)
 .map{ row -> tuple(row.sample, row.reads1, row.reads2)}
 .set { seedfile_ch }


/*
 * Run NinjaMap  16 128 fischbachlab/nf-ninjamap:20220726210531
 */
process ninjaMap {
    tag "$sample"
    container "fischbachlab/nf-ninjamap:latest"
    cpus { 16 * task.attempt }
    memory { 128.GB * task.attempt }

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${output_path}", mode:'copy'
    input:
    tuple val(sample), val(reads1), val(reads2) from seedfile_ch
    //file read1 from read1_ch
    //file read2 from read2_ch

    output:
    //path "*"

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
    ninjaMap_index_5.sh
    """
}

/*
process cal_depth {

  container "fischbachlab/nf-ninjamap:latest"
  cpus 4
  memory 8.GB

  input:
  file inbam from pe_bam_ch
  file abund from ninja_map_ch

  output:


  script:
  """
  s=${inbam%.bam}
  pileup.sh in=${inbam} out=${s}_coverage.txt overwrite=t 2>${s}_stats.txt
  echo -e "stain_name\tsingular_depth\tsingular_coverage" > ${s}_summary_depth.tsv
  for i in $(cut -f1 ${s}_coverage.txt | tail -n+2  | sed  's/_Node.*/
  /*      /' |  uniq )
  do
   grep $i ${s}_coverage.txt | awk '{t+=$2*$3; c+=$6; l+=$3}END{print $1"\t"t/l"\t"c/l}' >> ${s}_summary_depth.tsv
  done
  cut -f2,3 ${s}_summary_depth.tsv | paste ${abund} -  | awk '{print $1","$2","$3}' > tmp.ninjaMap.abundance.csv

  """

}
*/
