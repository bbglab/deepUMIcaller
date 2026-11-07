process COMPUTE_OPTIMAL_SEQUENCING {
    tag "compile"
    label 'family_metrics'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    path(duplex_metrics_all), stageAs: "duplex_metrics_all.tsv"
    path(curve_data_all), stageAs: "curve_data_all.tsv"
    path(duplex_metrics_on_target), stageAs: "duplex_metrics_on_target.tsv"
    path(curve_data_on_target), stageAs: "curve_data_on_target.tsv"


    output:
    path("*.pdf")              , optional: true ,  emit: pdf
    path("*.sample_data.tsv")  , optional: true ,  emit: sample_data
    path("*.family_curve.tsv") , optional: true ,  emit: curve_data
    path  "versions.yml"                        , topic: versions

    script:
    """
    metrics_compute_optimal.py \\
                --on_target_metrics ${duplex_metrics_on_target} \\
                --all_metrics ${duplex_metrics_all} \\
                --on_target_curve ${curve_data_on_target} \\
                --all_curve ${curve_data_all} \\
                --output-prefix compiled_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch compiled_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
