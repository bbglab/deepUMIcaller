process FGUMI_RETAGFROMCRAM {
    tag "$meta.id"
    label 'groupreads_io'

    conda "bioconda::fgumi bioconda::samtools=1.24"
    container 'fgumi:v0.8.0'

    input:
    tuple val(meta), path(cram), path(crai)
    path(fasta)

    output:
    tuple val(meta), path("*.retagged.bam"), emit: bam
    tuple val(meta), path("*.retagged.bam.bai"), emit: bai
    path "versions.yml"                      , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    samtools view -b \
        -T ${fasta} \
        -@ ${task.cpus} \
        -o ${prefix}.aligned.bam \
        ${cram}

    fgumi retag \
        --input ${prefix}.aligned.bam \
        --output ${prefix}.retagged.bam \
        rb,mb::pair::RX \
        rb::delete \
        mb::delete \
        --threads ${task.cpus} \
        $args

    samtools index -@ ${task.cpus} ${prefix}.retagged.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
        samtools: \$(samtools --version |& sed -n '1p' | sed 's/samtools //')
    END_VERSIONS
    """
}
