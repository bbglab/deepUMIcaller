process ALIGN_BAM {
    tag "$meta.id"
    label 'alignment_intensive'

    conda "bioconda::fgumi bioconda::bwa-mem3 bioconda::samtools"
    container 'wave.seqera.io/wt/4fd912ccf535/wave/build:fgumi_bwa-mem3_samtools--aefea4ff0576cc36'

    input:
    tuple val(meta), path(unmapped_bam)
    path index_dir

    output:
    tuple val(meta), path("*.mapped.bam")       , emit: bam
    tuple val(meta), path("*.mapped.bam.bai")   , emit: bai
    path "versions.yml"                         , topic: versions


    script:
    def bwa_args = task.ext.bwa_args ?: ''
    def fgumi_fastq_args = task.ext.fgumi_fastq_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def memory_gb = task.memory.toGiga() / 2
    def sort_cpus = task.cpus / 2
    """
    # The real path to the FASTA
    FASTA=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    mkdir temp_sort_directory

    fgumi fastq --input ${unmapped_bam} --threads ${task.cpus} ${fgumi_fastq_args} \\
        | bwa-mem3 mem ${bwa_args} -t $task.cpus -p -Y \$FASTA - \\

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
            --sort-threads ${sort_cpus} \
            --merge-threads ${task.cpus} \
            --tmp-dir temp_sort_directory/

    rm -rf temp_sort_directory

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem3: \$(bwa-mem3 version 2>&1 | head -1)
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
        bwa-mem3: \$(bwa-mem3 version 2>&1 | head -1)
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """

}
