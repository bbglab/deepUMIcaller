process UMICOLLISIONS {
    tag "$meta.id"
    label 'collect_duplex_metrics'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(position_group_sizes), path(umi_counts), path(duplex_umi_counts)

    output:
    tuple val(meta), path("*.umi_collisions.tsv")   , emit: tsv
    path "versions.yml"                             , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"

    // Handle single file or multiple files (list/collection), e.g. when split_by_chrom is enabled
    def position_group_files = [position_group_sizes].flatten().collect { file -> "--position-group-sizes ${file}" }.join(' ')
    def umi_counts_files = [umi_counts].flatten().collect { file -> "--umi-counts ${file}" }.join(' ')
    def duplex_umi_counts_files = [duplex_umi_counts].flatten().collect { file -> "--duplex-umi-counts ${file}" }.join(' ')

    def plot_arg = task.ext.plot ? "--plot" : ""

    """
    compute_umi_collisions.py \\
        --sample-name ${meta.id} \\
        ${position_group_files} \\
        ${umi_counts_files} \\
        ${duplex_umi_counts_files} \\
        --output-file ${prefix}.umi_collisions.tsv \\
        ${plot_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.umi_collisions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
