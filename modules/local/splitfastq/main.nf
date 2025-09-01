process SPLITFASTQ {
    tag "$meta.id"
    label 'process_low'
    label 'process_low_memory'

    conda "bioconda::seqkit=2.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0' :
        'biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("split_fastq/*_{1,2}.part_*.fastq.gz"), emit: split_fastqs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read1 = fastqs[0]
    def read2 = fastqs.size() > 1 ? fastqs[1] : ''
    def split_parts = task.ext.split_parts ?: params.splitfastq_parts ?: 20
    """
    mkdir -p split_fastq
    seqkit split2 -p ${split_parts} -O split_fastq $args -1 ${read1} ${read2 ? "-2 ${read2}" : ""}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version | sed 's/seqkit v//')
    END_VERSIONS
    """
}