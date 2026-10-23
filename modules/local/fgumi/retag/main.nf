process FGUMI_RETAGFROMCRAM {
    tag "$meta.id"
    label 'groupreads_io'

    conda "bioconda::fgumi=0.8.0 bioconda::samtools=1.24"
    container 'community.wave.seqera.io/library/fgumi_samtools:f4eabe7668db9d70'


    input:
    tuple val(meta), path(cram)
    path(fasta)

    output:
    tuple val(meta), path("*.retagged.bam")     , emit: bam
    tuple val(meta), path("*.retagged.bam.bai") , emit: bai
    path "versions.yml"                         , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -b \\
        -T ${fasta} \\
        -@ ${task.cpus} \\
        ${cram} | \\
        fgumi retag \\
            --input - \\
            --output - \\
            rb,mb::pair::RX \
            rb::delete \
            mb::delete \\
            --threads ${task.cpus} \
            ${args} | \\
            fgumi sort -i - -o ${prefix}.retagged.bam --threads ${task.cpus} --order coordinate --write-index true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
        samtools: \$(samtools --version |& sed -n '1p' | sed 's/samtools //')
    END_VERSIONS
    """
}