process QUERY_TABIX {
    tag "$meta.id"
    label 'cpu_single'
    label 'time_low'
    label 'memory_medium'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(pileup), path(pileuptabix), path(vcfderived)

    output:
    tuple val(meta), path("*.mutated_positions.tsv.gz") , emit: mutated_tsv
    path  "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # now we also need to recount the variants
    tabix ${pileup} -@ $task.cpus -R ${vcfderived} | bgzip > ${prefix}.mutated_positions.tsv.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mutated_positions.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}


// tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix (htslib) //; s/Copyright.*\$//')
