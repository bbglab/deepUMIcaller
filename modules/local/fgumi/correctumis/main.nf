process FGUMI_CORRECTUMIS {
    tag "$meta.id"
    label 'groupreads_io'

    conda "bioconda::fgumi"
    container  'quay.io/biocontainers/fgumi:0.1.3--h54198d6_0' 

    input:
    tuple val(meta), path(bam), path(umi_file)
    val(max_mismatches)
    val(min_distance)

    output:
    tuple val(meta), path("*.corrected.bam"), emit: bam
    tuple val(meta), path("*.rejected.bam"), emit: rejects
    tuple val(meta), path("*.correct-umis-metrics.txt"), emit: metrics
    path "versions.yml" , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fgumi correct \
        --input ${bam} \
        --output ${prefix}.corrected.bam \
        --rejects ${prefix}.rejected.bam \
        --metrics ${prefix}.correct-umis-metrics.txt \
        --max-mismatches ${max_mismatches} \
        --min-distance ${min_distance} \
        --umi-files ${umi_file} \
        --threads ${task.cpus} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
