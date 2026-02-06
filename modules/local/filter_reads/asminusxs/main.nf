process ASMINUSXS {
    tag "$meta.id"
    label 'bam_processing_heavy'
    
    conda "bioconda::pysam=0.23.3 conda-forge::click"
    container 'docker.io/bbglab/pysam-0.23.3:latest'

    input:
    tuple val(meta), path(bam), path (bam_index)

    output:
    tuple val(meta), path("*.AS-XS_*.bam")          , emit: bam
    tuple val(meta), path("*.discarded_AS-XS_*.bam"), emit: discarded_bam
    path  "versions.yml"                            , topic: versions


    script:
    def threshold = task.ext.threshold ?: "50"
    def prefix = task.ext.prefix ?: ".filtered.AS-XS_${threshold}"
    prefix = "${meta.id}${prefix}"
    def prefix_discard = task.ext.prefix_discard ?: ".discarded_AS-XS_${threshold}"
    prefix_discard = "${meta.id}${prefix_discard}"
    
    """
    as_minus_xs.py \
        --input_bam ${bam} \
        --output_bam ${prefix}.bam \
        --output_bam_discarded ${prefix_discard}.bam \
        --threshold ${threshold} \
        --tthreads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.23.3
        click: \$(python3 -c "import click; print(click.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def prefix_discard = task.ext.prefix_discard ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram
    touch ${prefix_discard}.bam
    touch ${prefix_discard}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pysam: 0.21.0
    END_VERSIONS
    """
}