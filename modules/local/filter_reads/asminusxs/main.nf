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
    tuple val(meta), path("*.discarded_AS-XS_*.bam"), emit: discarded_bam
    path  "versions.yml"                  , topic: versions


    script:
    def threshold = task.ext.threshold ?: "50"
    def prefix = task.ext.prefix ?: ".filtered.AS-XS_${threshold}"
    prefix = "${meta.id}${prefix}"
    def prefix_discard = task.ext.prefix_discard ?: ".discarded_AS-XS_${threshold}"
    prefix_discard = "${meta.id}${prefix_discard}"
    
    // TODO think of reimplementing with click
    """
    as_minus_xs.py ${bam} ${prefix}.bam ${prefix_discard}.bam ${threshold}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.21.0
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def prefix_discard = task.ext.prefix_discard ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram
    touch ${prefix_discard}.bam
    touch ${prefix_discard}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.21.0
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
