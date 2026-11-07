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
    def memory_gb = task.memory.toGiga().intdiv(2) + 1
    def sort_cpus = task.cpus.intdiv(2) + 1
    def reserve_memory_gb = task.memory.toGiga().intdiv(5) + 1
    """
    mkdir temp_sort_directory;
    fgumi sort \\
        --input $bam --output - \\
        --order queryname \\
        --max-memory ${memory_gb}GiB \
        --memory-per-thread false \
        --memory-reserve ${reserve_memory_gb}GiB \\
        --threads ${task.cpus} \
        --sort-threads ${sort_cpus} \
        --merge-threads ${task.cpus} \\
        --tmp-dir temp_sort_directory/ \\
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
