process CALLING_VARDICT {
    tag "$meta.id"
    label 'process_high_cpus'

    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"    

    input:
    tuple val(meta) , path(bam), path(bam_index), path (targets_file)
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
    total_lines=\$(wc -l < ${targets_file})
    # Calculate the number of lines per chunk
    lines_per_chunk=\$(( (total_lines + ${task.cpus} - 1) / ${task.cpus} ))
    # Split the file into chunks with the calculated number of lines
    split -l \${lines_per_chunk} ${targets_file} chunk_

    echo "Bed splitted in chunks. Running vardict-java..."

    # Process each chunk in parallel
    for chunk in chunk_*; do
        vardict-java -G ${fasta_dir}/${fasta} \
            -N ${prefix} -b ${bam} \
            -c 1 -S 2 -E 3 -g 4 \
            $args \
            -th 1 \
            \$chunk > \${chunk}.raw.tsv &
    done
    
    # Wait for all parallel processes to finish
    wait

    echo "Vardict finished. Concatenating..."

    # Concatenate all genome TSV chunks
    cat chunk_*.raw.tsv > ${prefix}.raw.tsv

    echo "Concatenated. teststrandbias running..."

    MAX_JOBS=4
    current_jobs=0

    for chunk in chunk_*.raw.tsv; do
        (
            cat \$chunk | teststrandbias.R | var2vcf_valid.pl \
                -N ${prefix} $filter_args \
                > \${chunk}.genome.vcf"
        ) &

        ((current_jobs++))
        if [[ \$current_jobs -ge \$MAX_JOBS ]]; then
            wait -n
            ((current_jobs--))
        fi
    done   
    
    # Wait for all parallel processes to finish
    wait

    echo "Done. Concatenating..."

    # Concatenate all genome VCF chunks
    # Extract the header from the first VCF chunk
    grep "^#" chunk_aa.raw.tsv.genome.vcf > ${prefix}.genome.vcf

    for chunk in chunk_*.genome.vcf; do
        grep -v "^#" \$chunk >> ${prefix}.genome.vcf
    done

    echo "Done. AWK filtering..."

    # Apply the AWK filter to create the final VCF
    awk '\$5!="."' ${prefix}.genome.vcf > ${prefix}.vcf

    echo "Done. Removing tmp files..."

    # Cleanup intermediate files (optional)
    rm chunk_*

    echo "Done. Gzip results..."

    gzip ${prefix}.raw.tsv;
    gzip ${prefix}.genome.vcf;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: 1.8.3
    END_VERSIONS
    """
}



