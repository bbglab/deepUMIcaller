process FGBIO_CLIPBAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fgbio=2.0.2 bioconda::bwa=0.7.17 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
            'https://depot.galaxyproject.org/singularity/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' : 
            'biocontainers/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' }"


    input:
    tuple val(meta), path(bam)
    path(fasta)


    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio ClipBam] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    """  
    samtools sort -n -u $bam \
        | fgbio \\
            -Xmx${mem_gb}g \\
            --tmp-dir=. \\
            ClipBam \\
            --input /dev/stdin/ \\
            --ref ${fasta} \\
            $args \\
            --output ${prefix}.clipped.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

