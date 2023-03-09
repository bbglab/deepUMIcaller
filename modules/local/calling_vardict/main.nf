process CALLING_VARDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0' }"    

    input:
    tuple val(meta), path(bam)
    path bam_index
    path targets_file
    path fasta
    path fasta_dir

    output:
    tuple val(meta), path("*.vcf")                , emit: vcf
    tuple val(meta), path("*.tsv"), optional: true, emit: tsv
    path  "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    AF_THR=0.0
    vardict-java -G ${fasta_dir}/${fasta} \\
        -f \$AF_THR \\
        -N ${prefix} -b ${bam} \\
        -c 1 -S 2 -E 3 -g 4 \\
        -r 1 -m 8 -P 0 \\
        -o 1 \\
        -th ${task.cpus} \\
        ${targets_file} > ${prefix}.raw.tsv

    cat ${prefix}.raw.tsv \\
        | teststrandbias.R \\
        | var2vcf_valid.pl \\
        -N ${prefix} -E -f \$AF_THR -p 0 -m 8 -v 2 > ${prefix}.vcf

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict: 1.8.3
    END_VERSIONS
    """
}

// vardict-java -G ${fasta_dir}/${fasta} \\
//         -f \$AF_THR \\
//         -N ${prefix} -b ${bam} \\
//         -c 1 -S 2 -E 3 -g 4 \\
//         -r 1 -m 8 -P 0 \\
//         -th ${task.cpus} \\
//         ${targets_file} \\
//         | teststrandbias.R \\
//         | var2vcf_valid.pl \\
//         -N ${prefix} -E -f \$AF_THR -p 0 -m 8 -v 2 > ${prefix}.vcf


    // vardict-java -G ${fasta_dir}/${fasta} \\
    //     -f \$AF_THR \\
    //     -N ${prefix} -b ${bam} \\
    //     -c 1 -S 2 -E 3 -g 4 \\
    //     -r 1 -m 8 -P 0 -D \\
    //     -o 1 \\
    //     -th ${task.cpus} \\
    //     ${targets_file} > ${prefix}.raw.vcf

    // cat ${prefix}.raw.vcf \\
    //     | teststrandbias.R \\
    //     | var2vcf_valid.pl \\
    //     -N ${prefix} -E -f \$AF_THR -p 0 -m 8 -v 2 > ${prefix}.vcf

