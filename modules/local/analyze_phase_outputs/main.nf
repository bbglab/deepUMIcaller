process ANALYZE_PHASE_OUTPUTS {
    tag "$meta.id"
    label 'postprocess_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(variant_chain_support_tsv)

    output:
    tuple val(meta), path("*.copy_number_by_chromosome.tsv")           , emit: copy_number_by_chromosome
    tuple val(meta), path("*.copy_number_by_window.tsv")               , emit: copy_number_by_window
    tuple val(meta), path("*.copy_number_calls.tsv")                   , emit: copy_number_calls
    tuple val(meta), path("*.independent_mutation_candidates.tsv")     , emit: independent_mutation_candidates
    tuple val(meta), path("*.analysis_summary.txt")                    , emit: analysis_summary
    path  "versions.yml"                                                , topic: versions

    script:
    def suffix = task.ext.suffix ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    prefix = suffix != '' ? "${prefix}.${suffix}" : prefix
    def output_prefix = "${prefix}.readjusted.mutated_reads"
    """
    analyze_phase_outputs.py \\
        --input ${variant_chain_support_tsv} \\
        --output-prefix ${output_prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def suffix = task.ext.suffix ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    prefix = suffix != '' ? "${prefix}.${suffix}" : prefix
    def output_prefix = "${prefix}.readjusted.mutated_reads"
    """
    touch ${output_prefix}.copy_number_by_chromosome.tsv
    touch ${output_prefix}.copy_number_by_window.tsv
    touch ${output_prefix}.copy_number_calls.tsv
    touch ${output_prefix}.independent_mutation_candidates.tsv
    touch ${output_prefix}.analysis_summary.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
