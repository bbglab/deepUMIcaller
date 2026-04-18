process FGUMI_CALLDUPLEXCONSENSUSREADS {
    tag "$meta.id"
    cache 'lenient'
    label 'consensus_calling'

    conda "bioconda::fgumi"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/fgumi:0.3.0--h4327870_0' :
        'quay.io/biocontainers/fgumi:0.3.0--h4327870_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , topic: versions


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}.consensus${prefix}"

    """
    fgumi duplex \
        --input $bam \
        --output ${prefix}.bam \
        --threads ${task.cpus} \
        --read-name-prefix ${meta.id} \
        --read-group-id ${meta.id} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
