process ALIGN_BAM {
    tag "$meta.id"
    label 'alignment_intensive'

    conda "bioconda::fgumi=0.8.0 bioconda::bwa-mem3=0.15.0"
    container 'community.wave.seqera.io/library/bwa-mem3_fgumi:6069fbf2f2bcec9f'

    input:
    tuple val(meta), path(unmapped_bam)
    path index_dir
    path ref_fasta

    output:
    tuple val(meta), path("*.mapped.bam")       , emit: bam
    tuple val(meta), path("*.mapped.bam.bai")   , emit: bai
    path "versions.yml"                         , topic: versions


    script:
    def bwa_args = task.ext.bwa_args ?: ''
    def fgumi_fastq_args = task.ext.fgumi_fastq_args ?: ''
    def fgumi_zipper_args = task.ext.fgumi_zipper_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def memory_gb = task.memory.toGiga().intdiv(3) + 1
    def sort_cpus = task.cpus.intdiv(2) + 1
    def reference_filename = ref_fasta.name
    def reserve_memory_gb = task.memory.toGiga().intdiv(5) + 1
    def in_n_out_task_cpus = 1
    """
    # The real path to the FASTA
    FASTA=`find -L ./ -name "${reference_filename}.amb" | sed 's/.amb//'`

    mkdir temp_sort_directory

    fgumi fastq --input ${unmapped_bam} \\
            --threads ${in_n_out_task_cpus} ${fgumi_fastq_args} \
            --max-memory ${reserve_memory_gb}GiB \
            --memory-per-thread false \\
        | bwa-mem3 mem ${bwa_args} -t $task.cpus -p -Y \$FASTA - \\
        | fgumi zipper \
            --input /dev/stdin \\
            --unmapped ${unmapped_bam} \
            --reference \$FASTA \
            --threads ${in_n_out_task_cpus} \
            ${fgumi_zipper_args} \\
        | fgumi sort --input /dev/stdin \
            --output ${prefix}.mapped.bam \
            --max-memory ${memory_gb}GiB \
            --memory-per-thread false \
            --order coordinate \
            --write-index true \
            --threads ${task.cpus} \
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
