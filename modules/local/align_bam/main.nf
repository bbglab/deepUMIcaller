process ALIGN_BAM {
    tag "$meta.id"
    label 'alignment_intensive'

    conda "bioconda::fgumi bioconda::bwa=0.7.17"
    container 'community.wave.seqera.io/library/bwa_fgumi_samtools:910c3ff2dc301fbf'

    input:
    tuple val(meta), path(unmapped_bam)
    path index_dir
    val sort

    output:
    tuple val(meta), path("*.mapped.bam"), emit: bam
    path "versions.yml"                  , topic: versions


    script:
    def bwa_args = task.ext.bwa_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    # The real path to the FASTA
    FASTA=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    fgumi fastq --input ${unmapped_bam} --threads ${task.cpus} \
        | bwa mem ${bwa_args} -t $task.cpus -p -Y \$FASTA - \
        | fgumi zipper \
            --input /dev/stdin \
            --unmapped ${unmapped_bam} \
            --reference \$FASTA \
            --threads ${task.cpus} \
        | fgumi sort --input /dev/stdin \\
            --output ${prefix}.mapped.bam --order template-coordinate

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
