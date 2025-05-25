process CREATEBED_FROM_TSV {
    tag "$meta.id"
    label 'process_medium'
    
    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"


    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tail -n +2 ${tsv} | awk -F'\\t' '{print \$1, \$2, \$2}' OFS='\\t' > ${prefix}.positions.tsv
    bedtools merge -i ${prefix}.positions.tsv \\
                    $args \\
                    > ${prefix}.sequenced.bed;

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

