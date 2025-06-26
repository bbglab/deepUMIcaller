process PATCH_DEPTH {
    tag "$meta.id"
    label 'process_single_medium_memory'
    
    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    tuple val(meta), path(pileup_mutations), path(vcf)

    output:
    tuple val(meta), path("*.readjusted.vcf")     , emit: patched_vcf
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def suffix = task.ext.suffix ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    prefix = suffix != '' ? "${prefix}.${suffix}" : prefix
    """
    recompute_depth.py \\
            ${pileup_mutations} \\
            ${vcf} \\
            ${prefix}.readjusted.vcf \\
            ${suffix} \\
            ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.readjusted.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

    // recompute_depth.py \\
    //         --mpileup_file ${pileup_mutations} \\
    //         --vcf_file ${vcf} \\
    //         --output ${prefix}.readjusted.vcf