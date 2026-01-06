/* short read QC */
process read_qc{

    tag "$sample"
    container  params.container
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }

    errorStrategy 'retry'
    maxRetries 1

    publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "Logs/bbduk.log.txt"
    publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "Stats/adapter_trimming_stats_per_ref.txt"

    input:
    tuple val(sample), path(R1), path(R2)

    output:
    tuple val(sample), path("qc_fastq/read1_trimmed.fastq.gz"), path("qc_fastq/read2_trimmed.fastq.gz"), emit: qc_fastq_ch
    path  "Logs/bbduk.log.txt"
    path  "Stats/adapter_trimming_stats_per_ref.txt"

    script:
    """
    export sampleRate="${params.sampleRate}"
    export trimQuality="${params.minQuality}"
    export minLength="${params.minLength}"
    export fastq1="${R1}"
    export fastq2="${R2}"
    export SAMPLE_NAME="${sample}"
    export autoSampling="${params.auto_sampling}"
    
    read_qc.sh 
    """
}


/*
*	Align each filtered sample to each strain
*/

process bowtie2_alignment {

    tag "$sample"
    errorStrategy 'retry' 
    maxRetries 1

    container params.container
    cpus 32
    memory 64.GB 
    
    publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy'

	//publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "bowtie2/*.{bam,bai,gz}"
    //publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "Logs/read_mapping.log.txt"
    //publishDir "${params.output_path}/${params.db}/${params.project}/${sample}/", mode:'copy', pattern: "Stats/read_accounting.csv"

    input:
    path "dir/*"   // entire bowtie2_index directory
	tuple val(sample), path(fq1), path(fq2), path(raw1), path(raw2)

    output:
    tuple val(sample), path("bowtie2/${sample}.sortedByCoord.bam"), path("bowtie2/${sample}.sortedByCoord.bam.bai"), optional: true, emit: bam_ch
    tuple val(sample), path("bowtie2/${sample}_unmapped_include_overlap_R1.fastq.gz"), path("bowtie2/${sample}_unmapped_include_overlap_R2.fastq.gz"), optional: true
    path "Logs/read_mapping.log.txt", optional: true
    path "Stats/read_accounting.csv", optional: true
   
    script:
    """
    export SAMPLE_NAME="${sample}"
    export coreNum="${params.coreNum}"
    export memPerCore="${params.memPerCore}"
    export REFDBNAME="${params.db_prefix}"

    bowtie2.sh dir $fq1 $fq2 $raw1
    """
}
