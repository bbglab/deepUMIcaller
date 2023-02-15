process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium_mem'

    conda (params.enable_conda ? "bioconda::fgbio=2.0.2 bioconda::samtools=1.16.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.0.2--hdfd78af_0' :
        'quay.io/biocontainers/fgbio:2.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path fasta
    val(min_reads)
    val(min_baseq)
    val(max_base_error_rate)

    output:
    tuple val(meta), path("*.filtered.bam")       , emit: bam
    path "versions.yml"                           , emit: versions

    script:
    def fgbio_args = task.ext.fgbio_args ?: ''
    def samtools_args = task.ext.samtools_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FilterConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    // sort = false
    // if (sort) {
    //     fgbio_zipper_bams_output = "/dev/stdout"
    //     fgbio_zipper_bams_compression = 0 // do not compress if samtools is consuming it
    //     extra_command = " | samtools sort "
    //     extra_command += samtools_sort_args
    //     extra_command += " --template-coordinate"
    //     extra_command += " --threads "+ task.cpus
    //     extra_command += " -o " + prefix + ".filtered.bam##idx##"+ prefix + ".filtered.bam.bai"
    //     extra_command += " --write-index"
    //     extra_command += samtools_args
    // } else {
    fgbio_zipper_bams_output = prefix + ".filtered.bam"
    fgbio_zipper_bams_compression = 1
    extra_command = ""
    // }
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=${fgbio_zipper_bams_compression} \\
        FilterConsensusReads \\
        --input $grouped_bam \\
        --ref ${fasta} \\
        --min-reads ${min_reads} \\
        --min-base-quality ${min_baseq} \\
        --max-base-error-rate ${max_base_error_rate} \\
        --output ${fgbio_zipper_bams_output} \\
        $fgbio_args;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}
