process CALLING_VARDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"    

    input:
    tuple val(meta) , path(bam), path(bam_index)
    tuple val(meta2), path (targets_file)
    path fasta
    path fasta_dir

    output:
    tuple val(meta), path("*.vcf")                   , emit: vcf
    tuple val(meta), path("*.vcf.gz"), optional: true, emit: genome_vcf
    tuple val(meta), path("*.tsv.gz"), optional: true, emit: tsv
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def filter_args = task.ext.filter_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    # Split the targets file into ${task.cpus} chunks
    split -n l/${task.cpus} -d --additional-suffix=.targets ${targets_file} chunk_
    
    # Process each chunk in parallel
    for chunk1 in chunk_*.targets; do
        vardict-java -G ${fasta_dir}/${fasta} \
            -N ${prefix} -b ${bam} \
            -c 1 -S 2 -E 3 -g 4 \
            $args \
            -th 1 \
            \$chunk1 > \${chunk1}.raw.tsv &
    done
    
    # Wait for all parallel processes to finish
    wait

    # Concatenate all genome TSV chunks
    cat chunk_*.raw.tsv > ${prefix}.raw.tsv

    for chunk2 in chunk_*.raw.tsv; do
        (
            cat \$chunk2 \
            | teststrandbias.R \
            | var2vcf_valid.pl \
                -N ${prefix} $filter_args \
            > \${chunk2}.genome.vcf
        ) &
    done
    
    # Wait for all parallel processes to finish
    wait

    # Concatenate all genome VCF chunks
    cat chunk_*.genome.vcf > ${prefix}.genome.vcf
    
    # Apply the AWK filter to create the final VCF
    awk '\$5!="."' ${prefix}.genome.vcf > ${prefix}.vcf
    
    # Cleanup intermediate files (optional)
    rm chunk_*.raw.tsv
    rm chunk_*.genome.vcf
    
    gzip ${prefix}.raw.tsv;
    gzip ${prefix}.genome.vcf;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: 1.8.3
    END_VERSIONS
    """
}



