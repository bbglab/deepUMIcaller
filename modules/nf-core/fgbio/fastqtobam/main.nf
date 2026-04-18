process FGUMI_FASTQTOBAM {
    tag "$meta.id"
    label 'fastq_processing'

    conda "bioconda::fgumi"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/fgumi:0.3.0--h4327870_0' :
        'quay.io/biocontainers/fgumi:0.3.0--h4327870_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam") , emit: bam , optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    path "versions.yml"            , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def sample_name = args.contains("--sample") ? "" : "--sample ${prefix}"
    def library_name = args.contains("--library") ? "" : "--library ${prefix}"
    def output = prefix =~ /\.(bam|cram)$/ ? prefix : "${prefix}.bam"
    def inputs = reads.collect { it.toString() }.join(' ')
    def read_structures = meta.read_structure ?: "+T +T"
    """

    fgumi extract \
        ${args} \
        --inputs ${inputs} \
        --output ${output} \
        --read-structures ${read_structures} \
        ${sample_name} \
        ${library_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
