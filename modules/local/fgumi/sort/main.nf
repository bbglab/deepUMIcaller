process FGUMI_SORT {
    tag "$meta.id"
    label 'bam_sort_compress'

    conda "bioconda::fgumi"
    label 'fgumi_tools' 

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    path  "versions.yml"          , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    if ("$bam".contains(".sorted.")) {
        prefix = prefix.replace(".sorted", ".resorted")
    }
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    def memory_gb = task.memory.toGiga().intdiv(2) + 1
    def sort_cpus = task.cpus.intdiv(2) + 1
    def reserve_memory_gb = task.memory.toGiga().intdiv(5) + 1
    """
    mkdir temp_sort_directory
    fgumi sort --input $bam \
            --output ${prefix}.bam \\
            --max-memory ${memory_gb}GiB \
            --memory-per-thread false \
            --memory-reserve ${reserve_memory_gb}GiB \\
            --sort-threads ${sort_cpus} \
            --merge-threads ${task.cpus} \
            --tmp-dir temp_sort_directory/ \\
            ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
