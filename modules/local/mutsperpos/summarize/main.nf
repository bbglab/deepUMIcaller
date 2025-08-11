process SUMMARIZE_MUTS_PER_POS {
    tag "cohort"

    label 'process_single'
    label 'time_medium'
    label 'cpu_low'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    path (csv_files)

    output:
    path 'ratios_per_sample.tsv'                , optional : true,  emit: ratios_table
    path 'mutation_ratios_summary.pdf'          , optional : true,  emit: summary_pdf
    path 'samples_passing_ratio_threshold.tsv'  , optional : true,  emit: ratios_table
    path  "versions.yml"                        , topic: versions

    script:
    def args = task.ext.args ?: '--initial 5'
    """
    summarize_muts_per_cycle.py ${args}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch mutation_ratios_summary.pdf
    touch ratios_per_sample.tsv
    touch samples_passing_ratio_threshold.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
