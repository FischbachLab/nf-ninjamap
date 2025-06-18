/*
 * Run NinjaMap  16 128 fischbachlab/nf-ninjamap:20220726210531
 */
process ninjaMap {
    tag "$sample"
    container  params.container
    cpus  32          // { 32 * task.attempt }
    memory 250.GB     //{ 250.GB * task.attempt }

    errorStrategy 'retry'  //{ task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    //publishDir "${params.output_path}/${params.db_prefix}/${params.project}", mode:'copy'

    input:
    tuple val(sample), val(reads1), val(reads2)

    output:
    path "tmp_*/Sync/ninjaMap/${sample}.ninjaMap.abundance.csv"   //, emit: abund_ch
    //path "tmp_*/Sync/ninjaMap/${sample}.ninjaMap.read_stats.csv", emit: stats_ch
    //path "tmp_*/Sync/Logs/${sample}_Contaminants_stats.csv", emit: contam_ch
    //path "tmp_*/Sync/Stats/read_accounting.csv", emit: read_ch

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
    export S3OUTPUTPATH="${params.output_path}/${params.db}/${params.project}/${sample}"
    export STRAIN_MAP_FILENAME="${params.db_prefix}.ninjaIndex.binmap.csv"
    export trimQuality="${params.minQuality}"
    export minLength="${params.minLength}"
    export debug="${params.debug}"
    export ref_db="${params.ref_db_path}"
    ninjaMap.sh
    """
}

/*
 * Aggregate ninjaMap results into figures and tables from the project dir
 */

 process aggregate_ninjaMap {

    tag  params.project
    container  params.container_R
    cpus { 2 * task.attempt }
    memory { 8.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.output_path}/${params.db}/${params.project}", mode:'copy'

    input:
      path "ninjamap_abundance_list/*"
      //path "ninjamap_read_stats._list/*" 
      //path "ninjamap_Contaminants_list/*" 
      //path "read_accounting.csv"

    output:
      path "aggregated_ninjamap_results/*"

    script:
    """
    aggregate_ninjamap_results.R  '${params.output_path}' ${params.db} ${params.project}
    """
 }