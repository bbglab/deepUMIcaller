process ALIGN_BAM {
    tag "$meta.id"
    label 'alignment_intensive'

    conda "bioconda::fgumi bioconda::bwa=0.7.17"
    container 'biocontainers/fgumi:0.1.3--h54198d6_0'

    input:
    tuple val(meta), path(unmapped_bam)
    path index_dir
    val sort

    output:
    tuple val(meta), path("*.mapped.bam"), emit: bam
    path "versions.yml"                  , topic: versions


    script:
    def samtools_sort_args = task.ext.samtools_sort_args ?: ''
    def bwa_args = task.ext.bwa_args ?: ''
    def fgumi_args = task.ext.fgumi_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    if (sort) {
        fgumi_zipper_bams_output = "/dev/stdout"
        extra_command = " | fgumi sort --input /dev/stdin --order template-coordinate --output ${prefix}.mapped.bam --threads ${task.cpus}"
    } else {
        fgumi_zipper_bams_output = prefix + ".mapped.bam"
        extra_command = ""
    }

    """
    # The real path to the FASTA
    FASTA=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    fgumi fastq --input ${unmapped_bam} --threads ${task.cpus} \
        | bwa mem ${bwa_args} -t $task.cpus -p -Y \$FASTA - \
        | fgumi zipper \
            --input /dev/stdin \
            --unmapped ${unmapped_bam} \
            --reference \$FASTA \
            --output ${fgumi_zipper_bams_output} \
            --threads ${task.cpus} \
            --tags-to-reverse Consensus \
            --tags-to-revcomp Consensus \
            ${fgumi_args} \
            ${extra_command}

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
