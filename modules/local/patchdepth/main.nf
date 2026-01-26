process PATCH_DEPTH {
    tag "$meta.id"
    label 'postprocess_memory'

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"

    input:
    tuple val(meta), path(pileup_mutations), path(vcf)

    output:
    tuple val(meta), path("*.readjusted.vcf")     , emit: patched_vcf
    path  "versions.yml"                          , topic: versions


    script:
    def suffix = task.ext.suffix ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    prefix = suffix != '' ? "${prefix}.${suffix}" : prefix
    suffix = suffix != '' ? "--suffix ${suffix}" : ''
    """
    recompute_depth.py \\
            --mpileup-file ${pileup_mutations} \\
            --vcf-file ${vcf} \\
            --output-filename ${prefix}.readjusted.vcf \\
            ${suffix}
            

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
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