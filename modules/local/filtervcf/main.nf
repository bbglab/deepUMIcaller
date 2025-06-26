process FILTERMUTATIONS {
    tag "$meta.id"
    label 'process_low'
    
    // TODO
    // update this in the nfcore format once the container is available in biocontainers and galaxy singularity
    conda "anaconda::seaborn=0.12.2"
    container "biocontainers/seaborn:0.12.2_cv1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
    //         'https://depot.galaxyproject.org/singularity/seaborn:0.12.2_cv1' : 
    //         'biocontainers/seaborn:0.12.2_cv1' }"


    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.filter_mutations.vcf"), emit: vcf
    tuple val(meta), path("*.png") , optional:true , emit: png
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vaf_filter = task.ext.vaf_filter ?: ''
    def filters = task.ext.filters ?: ''
    """
    filtervcf.py \\
                ${prefix} \\
                ${vcf} \\
                ${filters} \\
                ${vaf_filter} \\
                ${prefix}.filter_mutations.vcf

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

