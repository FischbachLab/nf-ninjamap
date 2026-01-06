/*######################################################################
#  HOST contamination check
######################################################################
*/
process ninjaMap_contam {

    tag "$SAMPLE_NAME"
    container  params.container
    cpus { 16 * task.attempt }
    memory { 32.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}", mode:'copy'

    input:
     // tuple val(SAMPLE_NAME), path(read1_trimmed), path(read2_trimmed), val(tReads), val(uReads)
     tuple val(SAMPLE_NAME), path(read1_trimmed), path(read2_trimmed), path(bam), path(bai)
 
    output:
      path "Logs/Host_Contaminants_stats.csv"

    when:
        params.check_human && params.check_mouse
    script:
    """
    export SAMPLE_NAME="${SAMPLE_NAME}"
    export read1="${read1_trimmed}"
    export read2="${read2_trimmed}"
    export ref="${params.ref_db_path}"
    export thres="${params.alignmentThres}"
    
    contam.sh $bam $bai
    """
}


/*##################################################################
Add Singular and Escrow coverage depth to the abundance table
#Notet that the sum of Singular and Escrow coverage coverage is often over 100%
################################################################
  set to 1 if enabling the debug mode 
  coverage must set to 1 if enabling the debug mode
  for generate singular&escrow bam files for each genome
*/
process ninjaMap_coverage_debug {

    tag "$SAMPLE_NAME"
    container  params.container_Bio
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}/ninjaMap", mode:'copy'

	when:
	  //params.coverage == '1' // NOT WORKING

    input:
      tuple val(SAMPLE_NAME), path(singular_bam), path(escrow_bam), path(abundance_full)
      path(REF)

    output:
      tuple val(SAMPLE_NAME), path("${SAMPLE_NAME}.ninjaMap.abundance.csv"), emit: abundance
      path("${SAMPLE_NAME}.singular.bed"), optional: true
      path("${SAMPLE_NAME}.escrow.bed"), optional: true
      path("${SAMPLE_NAME}.singular_sorted.bam"), optional: true
      path("${SAMPLE_NAME}.singular_sorted.bam.bai"), optional: true
      path("${SAMPLE_NAME}.escrow_sorted.bam"), optional: true
      path("${SAMPLE_NAME}.escrow_sorted.bam.bai"), optional: true
      path("debug"), optional: true

    script:
    """
    export SAMPLE_NAME="${SAMPLE_NAME}"
    export debug="${params.debug}"
    export coverage="${params.coverage}"
    export REF="${REF}"

    export in_abundance="$abundance_full"
    export singular_bam="$singular_bam"
    export escrow_bam="$escrow_bam"

    coverage_debug.sh
    """
}
/*
##################################################################
Compute Singular and Escrow coverage & depth in parallel,
then add them to the abundance table
################################################################
*/
process ninjaMap_coverage{

    tag "$SAMPLE_NAME"
    container  params.container_Bio
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}/ninjaMap", mode:'copy'

    input:
      tuple val(SAMPLE_NAME), path(singular_bam), path(escrow_bam), path(abundance_full)
      path(REF)

    output:
      tuple val(SAMPLE_NAME), path("${SAMPLE_NAME}.ninjaMap.abundance.csv"), emit: abundance
      path("${SAMPLE_NAME}.singular.bed"), optional: true
      path("${SAMPLE_NAME}.escrow.bed"), optional: true
      path("${SAMPLE_NAME}.singular_sorted.bam"), optional: true
      path("${SAMPLE_NAME}.singular_sorted.bam.bai"), optional: true
      path("${SAMPLE_NAME}.escrow_sorted.bam"), optional: true
      path("${SAMPLE_NAME}.escrow_sorted.bam.bai"), optional: true

    script:
    """
    export SAMPLE_NAME="${SAMPLE_NAME}"
    export coverage="${params.coverage}"
    export REF="${REF}"
    export in_abundance="$abundance_full"
    export singular_bam="$singular_bam"
    export escrow_bam="$escrow_bam"

    coverage_parallel.sh
    """
}
/*
##############################################################
Run the debug mode in parallel
coverage must set to 1 if enabling the debug mode
for generate singular&escrow bam files for each genome
################################################################  
*/
process ninjaMap_debug {

    tag "$SAMPLE_NAME"
    container  params.container_Bio
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}/ninjaMap", mode:'copy'

	when:
	  //params.debug == '1' // NOT WORKING

    input:
      tuple val(SAMPLE_NAME), path(singular_bam), path(escrow_bam), path(abundance_full)
      path(REF)

    output:
      path("debug"), optional: true

    script:
    """
    export SAMPLE_NAME="${SAMPLE_NAME}"
    export debug="${params.debug}"
    export REF="${REF}"
    export singular_bam="$singular_bam"
    export escrow_bam="$escrow_bam"

    debug_parallel.sh
    """
}



/*
Print a file on completion
*/
process ninjaMap_complete{

    tag  params.project
    container  params.container
    cpus {  2 * task.attempt }
    memory { 4.GB * task.attempt }

    errorStrategy 'retry'   // { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${SAMPLE_NAME}", mode:'copy',  pattern: "job.complete"

    input:
        tuple val(SAMPLE_NAME), path(abundance)

    output:
        path("job.complete")
        path("${SAMPLE_NAME}.txt"), emit: txt

    script:
    """
    date
  
    printf 'NinjaMap completed\n' > job.complete
    echo "Live long and prosper" >> job.complete
    echo "Done" > "${SAMPLE_NAME}.txt"
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
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}", mode:'copy'

    when: 
       params.plot

    input:
      path "ninjamap_complete_list/*"

    output:
      path "aggregated_ninjamap_results/*"

    script:
    """
    aggregate_ninjamap_results.R  '${params.output_path}' ${params.db} ${params.project}
    """
 }