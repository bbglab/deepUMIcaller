process FAMILYSIZEMETRICS {
    tag "$meta.id"
    errorStrategy 'ignore'
    label 'process_medium_low'

    // TODO
    // update this in the nfcore format once the container is available in biocontainers and galaxy singularity
    conda "anaconda::seaborn=0.12.2"
    container "biocontainers/seaborn:0.12.2_cv1"


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
    """
    family_size_plots.py \\
                ${sample_name} \\
                ${prefix}.duplex_seq_metrics.duplex_family_sizes.txt \\
                ${prefix}.family_sizes_plot_n_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.family_sizes_plot_n_stats.pdf \\ 
            ${prefix}.family_sizes_plot_n_stats.high.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

