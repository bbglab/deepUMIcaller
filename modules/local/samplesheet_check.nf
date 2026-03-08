process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'genomic_prep'
    
    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    val step

    output:
    path '*.csv'            , emit: csv
    path 'splitted'         , emit: splitted_input
    path "versions.yml"     , topic: versions

    script:
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv --log-level DEBUG --step $step

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
