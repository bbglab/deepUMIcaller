process FGUMI_CLIPBAM {
    tag "$meta.id"
    label 'bam_processing_heavy'

    conda "bioconda::fgumi"
    label 'fgumi_tools' 

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
    // FIXME : this should be switched to fgumi sort
    """
    fgumi sort \\
        --order queryname \
        --threads ${task.cpus} \
        --input $bam --output - \\
        | fgumi clip \\
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
