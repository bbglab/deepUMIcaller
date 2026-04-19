process FGUMI_GROUPREADSBYUMI {
    tag "$meta.id"
    cache 'lenient'
    label 'groupreads_io'

    conda "bioconda::fgumi"
    container  'quay.io/biocontainers/fgumi:0.1.3--h54198d6_0' 

    input:
    tuple val(meta), path(taggedbam)
    val(strategy)

    output:
    tuple val(meta), path("*_umi-grouped.bam")  , emit: bam
    tuple val(meta), path("*_umi_histogram.txt"), emit: histogram
    path "versions.yml"                         , topic: versions


    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def strategy_lc = strategy.toString().toLowerCase()

    """
    fgumi group \
        $args \
        --strategy ${strategy_lc} \
        --input $taggedbam \
        --output ${prefix}_umi-grouped.bam \
        --family-size-histogram ${prefix}_umi_histogram.txt \
        --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
