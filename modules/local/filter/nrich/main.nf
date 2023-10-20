// /workspace/projects/bladder_ts/scripts/add_filters_vcf

process FILTER_N_RICH {
    tag "$meta.id"
    label 'cpu_single'
    label 'time_low'
    label 'memory_medium'
    
    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
            'https://depot.galaxyproject.org/singularity/pandas:1.5.2' : 
            'biocontainers/pandas:1.5.2' }"


    input:
    tuple val(meta), path(vcf_file), path(ns_position_file), path(ns_position_index)

    output:
    tuple val(meta), path("*.filtered.vcf"), emit: filtered_vcf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO
    // add the depth limit as a parameter, also test how it works
    // also it would be better to set it dynamically as
    // a given quantile of depth? (25%?)
    """
    add_filter_nrich.py \\
            ${vcf_file} \\
            ${ns_position_file} \\
            ${prefix}.filtered.vcf \\
            n_rich \\
            1000
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

