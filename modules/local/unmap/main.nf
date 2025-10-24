process UNMAP_BAM {
    tag "$meta.id"
    label 'process_high'
    label 'process_medium_high_memory'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.noPG.bam") , emit: bam
    path  "versions.yml"                , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    if ("$bam" == "${prefix}.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools sort -@ ${task.cpus-1} -n -o ${prefix}.queryname.bam ${bam}
    samtools view -@ ${task.cpus-1} -H ${prefix}.queryname.bam | grep -v '^@PG' > header.tmp
    samtools reheader header.tmp ${prefix}.queryname.bam > ${prefix}.queryname.noPG.bam
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
    // samtools \\
    //     view \\
    //     --threads ${task.cpus-1} \\
    //     $args \\
    //     ${reference} \\
    //     ${readnames} \\
    //     -o ${prefix}.${file_type} \\
    //     $input \\
    //     $args2
