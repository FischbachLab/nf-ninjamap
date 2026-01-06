/*
 * Run NinjaMap  16 128 fischbachlab/nf-ninjamap:20220726210531
 */
process ninjaMap_abundance {
    tag "$sample"
    container  params.container
    cpus  32        //48 { 32 * task.attempt }
    memory 480.GB  // 250.GB     // 380 { 250.GB * task.attempt }

    errorStrategy 'retry'  //{ task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.output_path}/${params.db}/${params.project}/${sample}", mode:'copy'
    //publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "ninjaMap/*.ninjaMap.votes.csv.gz"
    //publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/debug/", mode:'copy', pattern: "*.csv"

    input:
    tuple val(sample), path(bam), path(bai)
    path (binmap)

    output:
    tuple val(sample), path("ninjaMap/${sample}.ninjaMap.abundance.csv"), optional: true, emit: abundance
    // The entire tuple must be optional or not optional
    tuple val(sample), path("ninjaMap/${sample}.singular.bam"), path("ninjaMap/${sample}.escrow.bam"), optional: true, emit: abundance_coverage
    path "ninjaMap/*.ninjaMap.votes.csv.gz", optional: true
    path "ninjaMap/*.csv", optional: true

    script:
    """
    export sampleRate="${params.sampleRate}"
    export coreNum="${params.coreNum}"
    export memPerCore="${params.memPerCore}"
    export coverage="${params.coverage}"
    export REFDBNAME="${params.db_prefix}"
    export STRAIN_MAP_FILENAME="${binmap}"
    export trimQuality="${params.minQuality}"
    export minLength="${params.minLength}"
    export debug="${params.debug}"
    export ref_db="${params.ref_db_path}"
    export singular_vote="${params.min_singular_vote}"
    export mask_bed="${params.mask}"
    export OUTPUT_PREFIX="${sample}"
    bash ninjaMap.sh $bam $bai
    """
}

/*
 * Run NinjaMap  16 128 fischbachlab/nf-ninjamap:20220726210531
 */
process ninjaMap_coverage_general {
    tag "$SAMPLE_NAME"
    container  params.container
    cpus  8        
    memory 16.GB  

    errorStrategy 'retry'  //{ task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2
    
    //publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}/debug/", mode:'copy'

    input:
    tuple val(SAMPLE_NAME), path(bam), path(bai)
    path (REF)

    output:
    tuple val(SAMPLE_NAME), path("${SAMPLE_NAME}_summary_coverage.tsv"), emit: coverage

    script:
    """
    coverage.sh $SAMPLE_NAME $bam $bai $REF
   
    """
}

process ninjaMap_merge_abundance_table{

    tag "$SAMPLE_NAME"
    container  params.container
    cpus  2       
    memory 4.GB  

    errorStrategy 'retry'  //{ task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}/ninjaMap", mode:'copy'

    input:
    // NO optional input file
    tuple val(SAMPLE_NAME), path(abundance), path(summary_coverage)

    output:
    tuple val(SAMPLE_NAME), path("${SAMPLE_NAME}.ninjaMap.abundance.csv"), emit: abundance
    //path ("merged_output.csv")

    script:
    """
    merge_abundance_table.py $SAMPLE_NAME $abundance $summary_coverage
    """
}