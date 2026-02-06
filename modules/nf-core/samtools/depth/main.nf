process SAMTOOLS_DEPTH {
    tag "$meta.id"
    label 'coverage_compute'
    
    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        -o ${prefix}.tsv \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
