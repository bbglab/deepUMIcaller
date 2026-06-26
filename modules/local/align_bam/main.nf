process ALIGN_BAM {
    tag "$meta.id"
    label 'alignment_intensive'

    conda "bioconda::fgumi bioconda::bwa=0.7.17"
    container 'community.wave.seqera.io/library/bwa_fgumi_samtools:86e1d6ef7afef498'

    input:
    tuple val(meta), path(unmapped_bam)
    path index_dir

    output:
    tuple val(meta), path("*.mapped.bam")       , emit: bam
    tuple val(meta), path("*.mapped.bam.bai")   , emit: bai
    path "versions.yml"                         , topic: versions


    script:
    def bwa_args = task.ext.bwa_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def memory_gb = task.memory.toGiga() / 2
    """
    # The real path to the FASTA
    FASTA=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    mkdir temp_sort_directory

    fgumi fastq --input ${unmapped_bam} --threads ${task.cpus} \\
        | bwa mem ${bwa_args} -t $task.cpus -p -Y \$FASTA - \\
        | fgumi zipper \
            --input /dev/stdin \
            --unmapped ${unmapped_bam} \
            --reference \$FASTA \
            --threads ${task.cpus} \\
        | fgumi sort --input /dev/stdin \
            --output ${prefix}.mapped.bam \
            --max-memory ${memory_gb}GiB \
            --memory-per-thread false \
            --memory-reserve 2GiB \
            --order coordinate \
            --write-index true \
            --threads ${task.cpus} \
            --tmp-dir temp_sort_directory/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """

    stub:

    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.mapped.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """

}
