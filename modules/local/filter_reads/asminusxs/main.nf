process ASMINUSXS {
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_memory'

    conda "bioconda::pysam-0.21.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam-0.21.0--py38h15b938a_1' :
        'biocontainers/pysam-0.21.0--py38h15b938a_1' }"

    input:
    tuple val(meta), path(bam), path (bam_index)

    output:
    tuple val(meta), path("*.AS-XS_*.bam"), emit: bam
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def threshold = task.ext.threshold ?: "50"
    def prefix = task.ext.prefix ?: "${meta.id}.filtered.AS-XS_${threshold}"
    
    """
    as_minus_xs.py ${bam} ${prefix}.bam ${threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
    // samtools \\
    //     view \\
    //     --threads ${task.cpus-1} \\
    //     $args \\
    //     ${reference} \\
    //     ${readnames} \\
    //     -o ${prefix}.${file_type} \\
    //     $input \\
    //     $args2
