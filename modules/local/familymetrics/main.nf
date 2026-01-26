process FAMILYSIZEMETRICS {
    tag "$meta.id"
    label 'family_metrics'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(duplex_metrics)

    output:
    tuple val(meta), path("*.pdf")              , emit: pdf
    tuple val(meta), path("*.sample_data.tsv")  , emit: sample_data
    tuple val(meta), path("*.family_curve.tsv") , emit: curve_data
    path  "versions.yml"                        , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    def sample_name = "${meta.id}"
    prefix = "${meta.id}${prefix}"
    def confidence = task.ext.confidence ? "--confidence-level '${task.ext.confidence}' " : ""
    
    // Handle single file or multiple files (list/collection)
    def input_files = [duplex_metrics].flatten().collect { file -> "--input-file ${file}" }.join(' ')    
    """
    family_size_plots.py \\
                --sample-name ${sample_name} \\
                ${input_files} \\
                --output-file ${prefix}.family_sizes_plot_n_stats \\
                ${confidence}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.family_sizes_plot_n_stats.pdf \ 
            ${prefix}.family_sizes_plot_n_stats.high.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
