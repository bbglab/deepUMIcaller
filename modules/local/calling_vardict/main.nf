process CALLING_VARDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"    

    input:
    tuple val(meta), path(bam), path(bam_index)
    path targets_file
    path fasta
    path fasta_dir

    output:
    tuple val(meta), path("*.vcf")                   , emit: vcf
    tuple val(meta), path("*.vcf.gz"), optional: true, emit: genome_vcf
    tuple val(meta), path("*.tsv")   , optional: true, emit: tsv
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def filter_args = task.ext.filter_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    vardict-java -G ${fasta_dir}/${fasta} \\
        -N ${prefix} -b ${bam} \\
        -c 1 -S 2 -E 3 -g 4 \\
        $args \\
        -th ${task.cpus} \\
        ${targets_file} > ${prefix}.raw.tsv

    cat ${prefix}.raw.tsv \\
        | teststrandbias.R \\
        | var2vcf_valid.pl \\
        -N ${prefix} $filter_args \\
        | gzip > ${prefix}.genome.vcf.gz
    
    zcat ${prefix}.genome.vcf.gz | awk '\$5!="."' > ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: 1.8.3
    END_VERSIONS
    """
}



