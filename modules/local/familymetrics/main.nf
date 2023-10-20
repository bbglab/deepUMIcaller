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
    tuple val(meta), path(groupby_metrics), path(duplex_metrics)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path  "versions.yml"          , emit: versions
    stdout                          emit: log


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    family_size_plots.py \\
                ${prefix} \\
                ${groupby_metrics} \\
                ${prefix}.duplex_seq_metrics.duplex_family_sizes.txt \\
                ${prefix}.family_sizes_plot_n_stats.pdf \\
                ${prefix}.family_sizes_plot_n_stats.high.pdf

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

