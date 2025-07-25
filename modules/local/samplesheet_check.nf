process SAMPLESHEET_CHECK {
    tag "$samplesheet"

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path samplesheet
    val step

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", topic: versions

    script: // This script is bundled with the pipeline, in nf-core/fgcons/bin/
    """
    check_samplesheet.py $samplesheet samplesheet.valid.csv --step $step

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
