process BEDTOOLS_COVERAGE {
    tag "$meta.id"
    label 'process_medium'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.1--hf5e1c6e_0' :
        'biocontainers/bedtools:2.31.1--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(input_A), path(input_B)
    path genome_file

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def reference = genome_file ? "-g ${genome_file} -sorted" : ""
    """
    bedtools \\
        coverage \\
        $args \\
        $reference \\
        -a $input_A \\
        -b $input_B \\
        > ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(echo \$(bedtools --version 2>&1) | sed 's/^.*bedtools v//' ))
    END_VERSIONS
    """
}