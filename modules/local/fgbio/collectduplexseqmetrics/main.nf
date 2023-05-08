process FGBIO_COLLECTDUPLEXSEQMETRICS {
    tag "$meta.id"
    label 'process_low_multicpu'

    conda "bioconda::fgbio=2.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0' :
        'biocontainers/fgbio:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(grouped_bam)

    output:
    tuple val(meta), path("*duplex_seq_metrics*.txt")                , emit: metrics
    tuple val(meta), path("*.pdf")                   , optional: true, emit: report
    path "versions.yml"                                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio CollectDuplexSeqMetrics] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --async-io=true \\
        --compression=1 \\
        CollectDuplexSeqMetrics \\
        --input $grouped_bam \\
        --output ${prefix}.duplex_seq_metrics \\
        $args;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    touch $prefix.duplex_seq_metrics.duplex_qc.pdf
    touch $prefix.duplex_seq_metrics.duplex_umi_counts.txt
    touch $prefix.duplex_seq_metrics.umi_counts.txt
    touch $prefix.duplex_seq_metrics.duplex_yield_metrics.txt
    touch $prefix.duplex_seq_metrics.duplex_family_sizes.txt
    touch $prefix.duplex_seq_metrics.family_sizes.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

}
