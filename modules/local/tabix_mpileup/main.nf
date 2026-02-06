process QUERY_TABIX {
    tag "$meta.id"
    label 'postprocess_compute'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(pileup), path(pileuptabix), path(vcfderived)

    output:
    tuple val(meta), path("*.mutated_positions.tsv.gz") , emit: mutated_tsv
    path  "versions.yml"                                , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    # now we also need to recount the variants
    tabix ${pileup} -R ${vcfderived} | bgzip > ${prefix}.mutated_positions.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.mutated_positions.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}