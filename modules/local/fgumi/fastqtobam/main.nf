process FGUMI_FASTQTOBAM {
    tag "$meta.id"
    label 'fastq_processing'

    conda "bioconda::fgumi"
    container  'quay.io/biocontainers/fgumi:0.1.3--h54198d6_0' 

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("*.unmapped.bam"), emit: bam
    path "versions.yml"                    , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def read_structure = "${meta.read_structure}"
    def inputs = fastqs.collect { it.toString() }.join(' ')
    """

    fgumi extract \
        --inputs ${inputs} \
        --output "${prefix}.unmapped.bam" \
        --read-structures ${read_structure} \
        --sample ${meta.sample} \
        --library ${meta.sample} \
        --threads ${task.cpus} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
