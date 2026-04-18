process FGUMI_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'consensus_filter'

    conda "bioconda::fgumi"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/fgumi:0.3.0--h4327870_0' :
        'quay.io/biocontainers/fgumi:0.3.0--h4327870_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path fasta

    output:
    tuple val(meta), path("*.filtered.bam")  , emit: bam
    path "versions.yml"                      , topic: versions

    script:
    def fgumi_args = task.ext.fgumi_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"

    """
    fgumi filter \
        --input $grouped_bam \
        --ref ${fasta} \
        --output ${prefix}.filtered.bam \
        --threads ${task.cpus} \
        $fgumi_args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
