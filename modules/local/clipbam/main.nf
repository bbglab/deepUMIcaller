process FGUMI_CLIPBAM {
    tag "$meta.id"
    label 'bam_processing_heavy'

    conda "bioconda::fgumi bioconda::samtools=1.16.1"
    container 'biocontainers/fgumi:0.1.3--h54198d6_0'


    input:
    tuple val(meta), path(bam)
    path(fasta)


    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def manual_clipping = task.ext.extra_clipping ?: ""
    """
    samtools sort -n -@ ${task.cpus} -u $bam \
        | fgumi clip \
            --input /dev/stdin \
            --reference ${fasta} \
            ${manual_clipping} \
            $args \
            --output ${prefix}.clipped.bam \
            --threads ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
