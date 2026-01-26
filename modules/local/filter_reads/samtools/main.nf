process SAMTOOLS_FILTER {
    tag "$meta.id"
    label 'consensus_filter'

    conda "bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.0x2.bam") , emit: bam ,    optional: true
    tuple val(meta), path("*.cram"), emit: cram,    optional: true
    tuple val(meta), path("*.sam") , emit: sam ,    optional: true
    tuple val(meta), path("*.bai") , emit: bai ,    optional: true
    tuple val(meta), path("*.csi") , emit: csi ,    optional: true
    tuple val(meta), path("*.crai"), emit: crai,    optional: true
    path  "versions.yml"           , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    //def memory_per_cpu = Math.floor(task.memory.toGiga() / task.cpus) as Integer
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    bam.getExtension()
    if ("$bam" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    // TODO fix this into a more formal samtools view module that we call a particular name and provide the arguments in the modules.config
    """
    samtools sort \\
            --threads ${task.cpus-1} \\
            -n \\
            ${bam} \\
            | samtools fixmate \\
                    --threads ${task.cpus-1} \\
                    -r - ${meta.id}.tmp.bam

    samtools view --threads ${task.cpus-1} ${meta.id}.tmp.bam ${args} > ${prefix}.bam
    
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
