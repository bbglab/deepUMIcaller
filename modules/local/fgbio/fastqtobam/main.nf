process FGBIO_FASTQTOBAM {
    tag "$meta.id"
    label 'process_low'
    label 'process_low_memory'

    conda "bioconda::fgbio=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0' :
        'biocontainers/fgbio:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastqs)

    output:
    tuple val(meta), path("*.unmapped.bam"), emit: bam
    path "versions.yml"                    , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def mem_gb = 1
    def read_structure = "${meta.read_structure}"
    if (!task.memory) {
        log.info '[fgbio FastqToBam] Available memory not known - defaulting to 1GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    """

    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        FastqToBam \\
        --input ${fastqs} \\
        --output "${prefix}.unmapped.bam" \\
        --read-structures ${read_structure} \\
        --sample ${meta.id} \\
        --library ${meta.id} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
