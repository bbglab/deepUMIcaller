process FGUMI_COLLECTDUPLEXSEQMETRICS {
    tag "$meta.id"
    label 'collect_duplex_metrics'

    conda "bioconda::fgumi"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/fgumi:0.3.0--h4327870_0' :
        'quay.io/biocontainers/fgumi:0.3.0--h4327870_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path (intervals_file)

    output:
    tuple val(meta), path("*duplex_family_sizes.txt")                , emit: family_sizes
    tuple val(meta), path("*duplex_seq_metrics*.txt")                , emit: metrics
    tuple val(meta), path("*.pdf")                   , optional: true, emit: report
    path "versions.yml"                                              , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def coords_file = intervals_file ? "--intervals ${intervals_file}" : ""

    """
    fgumi duplex-metrics \
        --input $grouped_bam \
        --output ${prefix}.duplex_seq_metrics \
        ${coords_file} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    touch $prefix.duplex_seq_metrics.duplex_qc.pdf
    touch $prefix.duplex_seq_metrics.duplex_umi_counts.txt
    touch $prefix.duplex_seq_metrics.umi_counts.txt
    touch $prefix.duplex_seq_metrics.duplex_yield_metrics.txt
    touch $prefix.duplex_seq_metrics.duplex_family_sizes.txt
    touch $prefix.duplex_seq_metrics.family_sizes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """

}
