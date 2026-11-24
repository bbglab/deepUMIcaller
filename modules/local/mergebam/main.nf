process MERGEBAM {
    tag "$meta.id"
    label 'genomic_prep'
    
    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.merged.bam"), path("*.merged.bam.bai"), emit: bam_bai  // Changed from bam_and_index to bam_bai
    path "versions.yml", emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Create sorted list of input BAMs
    ls *.bam | sort > bam_list.txt
    n_bams=\$(wc -l < bam_list.txt)

    if [ "\$n_bams" -eq 1 ]; then
        # Single input (already coordinate-sorted): nothing to merge, so avoid a
        # full decompress/recompress pass and just link it and index it.
        ln -s "\$(readlink -f \$(cat bam_list.txt))" ${prefix}.merged.bam
        samtools index -@ $task.cpus ${prefix}.merged.bam
    else
        # Merge using the file list; the index is written in the same pass
        samtools merge $args -@ $task.cpus --write-index -b bam_list.txt \\
            -o ${prefix}.merged.bam##idx##${prefix}.merged.bam.bai
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.merged.bam
    touch ${prefix}.merged.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}