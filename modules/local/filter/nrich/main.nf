// /workspace/projects/bladder_ts/scripts/add_filters_vcf

process FILTER_N_RICH {
    tag "$meta.id"
    label 'process_memory_intensive'
    
    conda "conda-forge::pandas=2.3.2 conda-forge::click"
    container 'docker.io/bbglab/pysam-0.23.3:latest'


    input:
    tuple val(meta), path(vcf_file), path(ns_position_file), path(ns_position_index)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: filtered_vcf
    path  "versions.yml"                   , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def filter_name = task.ext.filter_name ?: "n_rich"
    def minimum_depth = task.ext.minimum_depth ?: "25"

    // TODO add this in a proper documentation
    // right now the value `minimum_depth` sets the minimum depth that a position needs to have
    // so that it can be used for computing the median coverage of the sample
    // once the median coverage of the sample is computed we only use positions 
    // with a coverage above (median * 0.75) to then compute the proportion of Ns
    // and with the proportion of Ns of all those positions we estimate a
    // distribution from which we compute the mean + 2 std as the threshold for the max. proportion of Ns

    // TODO think if we want to reimplement this with click
    """
    add_filter_nrich.py \\
            ${vcf_file} \\
            ${ns_position_file} \\
            ${prefix}.filtered.vcf \\
            ${filter_name} \\
            ${minimum_depth}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

