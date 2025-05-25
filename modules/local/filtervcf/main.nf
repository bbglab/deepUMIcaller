process FILTERMUTATIONS {
    tag "$meta.id"
    label 'process_medium'
    
    // TODO
    // update this in the nfcore format once the container is available in biocontainers and galaxy singularity
    conda "anaconda::seaborn=0.12.2"
    container "biocontainers/seaborn:0.12.2_cv1"


    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filter_mutations.vcf")                             , emit: vcf
    tuple val(meta), path("*.filter_mutations.purine.vcf")      , optional:true , emit: pur_vcf
    tuple val(meta), path("*.filter_mutations.pyrimidine.vcf")  , optional:true , emit: pyr_vcf
    tuple val(meta), path("*.png")                              , optional:true , emit: png
    path "versions.yml"                                                         , topic: versions


    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vaf_filter = task.ext.vaf_filter ?: ''
    def filters = task.ext.filters ?: ''
    def splitting = task.ext.splitting ?: ''
    """
    filtervcf.py \\
                ${prefix} \\
                ${vcf} \\
                ${filters} \\
                ${vaf_filter} \\
                ${prefix} \\
                ${splitting}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filter_mutations.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

