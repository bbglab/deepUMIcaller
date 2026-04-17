process PHASE_MUTATIONS {
    tag "$meta.id"
    label 'postprocess_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(mutated_reads_tsv)

    output:
    tuple val(meta), path("*.read_chains.tsv")                , emit: read_chains
    tuple val(meta), path("*.variant_chain_support.tsv")      , emit: variant_chain_support
    tuple val(meta), path("*.outliers.tsv")                   , emit: outliers
    tuple val(meta), path("*.pdf")          , optional: true  , emit: chain_proportions_plot
    path  "versions.yml"                                      , topic: versions

    script:
    def suffix = task.ext.suffix ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    prefix = suffix != '' ? "${prefix}.${suffix}" : prefix
    def output_prefix = "${prefix}.readjusted.mutated_reads"
    """
    phase_mutations.py \\
        --input ${mutated_reads_tsv} \\
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
    touch ${output_prefix}.read_chains.tsv
    touch ${output_prefix}.variant_chain_support.tsv
    touch ${output_prefix}.outliers.tsv
    touch ${output_prefix}.chain_proportions.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
