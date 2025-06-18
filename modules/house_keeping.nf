import groovy.json.JsonOutput
/*
  Save all parameters to a jason file
  https://github.com/nextflow-io/nextflow/discussions/2892
  //echo ${params} > ninjamap-parameters.txt
*/
process printParams {
    tag  params.project
    container  params.container
    publishDir "${params.output_path}/${params.db}/${params.project}/pipeline_info"

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
    tag  params.project
    container  params.container
    errorStrategy 'ignore'
    publishDir "${params.output_path}/${params.db}/${params.project}/pipeline_info"
    //saveAs: {filename ->
    //    if (filename.indexOf(".csv") > 0) filename
    //    else null
    //}

    output:
    path 'software_versions_ninjamap.txt'
    //path 'software_versions_ninja.txt'
    //path "software_versions.csv"
    //path "*.txt"

    script:
    // Get all tools to print their version number here
    // scrape_software_versions.py &> software_versions_ninjamap.txt
     //echo $workflow.manifest.version > v_pipeline.txt
    //echo $workflow.nextflow.version > v_nextflow.txt
    //python --version > v_python.txt
    //printf "pipeline_hash: %s\n" ${workflow.scriptId} >> software_versions_ninjamap.txt
    """
    printf "nextflow_version: %s\n" ${workflow.nextflow.version} > software_versions_ninjamap.txt
    printf "pipeline_version: %s\n" ${workflow.manifest.version} >> software_versions_ninjamap.txt
    printf "samtools_version: %s\n" \$(samtools --version | head -1 | awk '{print \$NF}') >> software_versions_ninjamap.txt
    printf "bedtools_version: %s\n" \$(bedtools --version | head -1 | awk -F " v" '{print \$2}') >> software_versions_ninjamap.txt
    printf "bowtie2_version: %s\n" \$(bowtie2 --version | grep -a bowtie2-align-s | awk '{print \$NF}') >> software_versions_ninjamap.txt
    printf "python_version: %s\n" \$(python --version | awk '{print \$NF}') >> software_versions_ninjamap.txt
    printf "biopython_version: %s\n" \$(python -c "import Bio; print(Bio.__version__)") >> software_versions_ninjamap.txt
    """
}

