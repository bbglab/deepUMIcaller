process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'consensus_filter'
    
    conda "bioconda::fgbio=2.1.0 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0' :
        'biocontainers/fgbio:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path fasta

    output:
    tuple val(meta), path("*.filtered.bam")  , emit: bam
    path "versions.yml"                      , topic: versions

    script:
    def fgbio_args = task.ext.fgbio_args ?: ''
    def Ns_per_read = (params.maxN_per_read + params.left_clip + params.right_clip) / params.read_length_minus_tag
    def max_no_call_fraction = Ns_per_read ? "--max-no-call-fraction ${Ns_per_read}" : "--max-no-call-fraction 0.2"
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FilterConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    fgbio_zipper_bams_output = prefix + ".filtered.bam"
    fgbio_zipper_bams_compression = 1
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=${fgbio_zipper_bams_compression} \\
        FilterConsensusReads \\
        --input $grouped_bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        ${max_no_call_fraction} \\
        ${fgbio_args};

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}