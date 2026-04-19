process FGUMI_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'consensus_filter'

    conda "bioconda::fgumi"
    container  'quay.io/biocontainers/fgumi:0.1.3--h54198d6_0' 

    input:
    tuple val(meta), path(bam)
    path fasta

    output:
    tuple val(meta), path("${prefix}.bam"), emit: bam
    path "versions.yml"                   , topic: versions


    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}.consensus_filtered${prefix}"

    """
    fgumi filter \
        --input $bam \
        --output ${prefix}.bam \
        --ref ${fasta} \
        --threads ${task.cpus} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
