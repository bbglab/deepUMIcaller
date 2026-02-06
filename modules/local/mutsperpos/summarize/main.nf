process SUMMARIZE_MUTS_PER_POS {
    tag "cohort"
    label 'cohort_analysis'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    path (csv_files)

    output:
    path 'ratios_per_sample_mqc.tsv'                , optional : true,  emit: ratios_table
    path 'mutation_ratios_summary_mqc.pdf'          , optional : true,  emit: summary_pdf
    path 'samples_passing_ratio_threshold_mqc.tsv'  , optional : true,  emit: failing_samples
    path  "versions.yml"                            , topic: versions

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
    touch mutation_ratios_summary_mqc.pdf
    touch ratios_per_sample_mqc.tsv
    touch samples_passing_ratio_threshold_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
