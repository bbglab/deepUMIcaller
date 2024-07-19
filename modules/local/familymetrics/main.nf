process FAMILYSIZEMETRICS {
    tag "$meta.id"
    label 'process_medium'
    
    // TODO
    // update this in the nfcore format once the container is available in biocontainers and galaxy singularity
    conda "anaconda::seaborn=0.12.2"
    container "biocontainers/seaborn:0.12.2_cv1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    //         'https://depot.galaxyproject.org/singularity/seaborn:0.12.2_cv1' : 
    //         'biocontainers/seaborn:0.12.2_cv1' }"


    input:
    tuple val(meta), path(duplex_metrics)

    output:
    tuple val(meta), path("*.pdf")              , emit: pdf
    tuple val(meta), path("*.sample_data.tsv")  , emit: sample_data
    tuple val(meta), path("*.family_curve.tsv") , emit: curve_data
    path  "versions.yml"                        , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def high_confidence = task.ext.high_confidence ?: "2 1 1"
    def med_confidence = task.ext.med_confidence ?: "2 1 1"
    def low_confidence = task.ext.low_confidence ?: "2 1 1"
    """
    family_size_plots.py \\
                ${prefix} \\
                ${prefix}.duplex_seq_metrics.duplex_family_sizes.txt \\
                ${prefix}.family_sizes_plot_n_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.family_sizes_plot_n_stats.pdf \\ 
            ${prefix}.family_sizes_plot_n_stats.high.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

