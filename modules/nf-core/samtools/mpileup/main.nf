process SAMTOOLS_MPILEUP {
    tag "$meta.id"
    label 'time_low'
    label 'memory_medium'

    conda "bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"
    input:
    tuple val(meta), path(bam), path(index), path(intervals)
    path fasta

    output:
    tuple val(meta), path("*.mpileup.gz"), path("*.mpileup.gz.tbi") , emit: mpileup
    path  "versions.yml"                                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = intervals ? "--positions ${intervals}" : ""
    """
    samtools mpileup \\
        --fasta-ref ${fasta} \\
        --output ${prefix}.mpileup \\
        ${args} \\
        ${intervals} \\
        ${bam}
    bgzip -@${task.cpus} ${prefix}.mpileup
    tabix -s 1 -b 2 -e 2 ${prefix}.mpileup.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
// samtools mpileup 
//         --no-BAQ
//         --max-depth 0
//         --fasta-ref $reference_fasta 
//         --positions merged_bedfile
//         --min-BQ 2
//         --no-output-ends 
//         $file
//         | \
//         bgzip > $samp.output_pileup.allowNs.notoverlapping.gz;