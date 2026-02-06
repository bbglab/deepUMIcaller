process MERGE_VARDICT_RESULTS {
    tag "$meta.id"
    label 'variant_calling'
    
    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf_chunks), path(raw_chunks)

    output:
    tuple val(meta), path("*.vcf")    , emit: vcf
    tuple val(meta), path("*.vcf.gz") , emit: genome_vcf
    tuple val(meta), path("*.tsv.gz"), optional: true, emit: tsv
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    echo "Merging VCF chunks..."

    # Sort chunks to ensure consistent ordering (chunk_aa, chunk_ab, etc.)
    vcf_files=(\$(ls *.genome.vcf | sort))

    # Extract the header from the first VCF chunk
    grep "^#" "\${vcf_files[0]}" > ${prefix}.genome.vcf

    # Concatenate all genome VCF chunks (excluding headers)
    for vcf_file in "\${vcf_files[@]}"; do
        grep -v "^#" "\$vcf_file" >> ${prefix}.genome.vcf
    done

    echo "Done concatenating VCFs. Applying AWK filter..."

    # Apply the AWK filter to create the final VCF
    awk '\$5!="."' ${prefix}.genome.vcf > ${prefix}.vcf

    echo "Done. Compressing genome VCF..."

    gzip ${prefix}.genome.vcf

    echo "Merging RAW TSV chunks..."

    # Merge raw TSV chunks (order by chunk name for determinism)
    raw_files=( \$(ls *.raw.tsv | sort) )
    if [ \${#raw_files[@]} -gt 0 ]; then
        cat "\${raw_files[@]}" > ${prefix}.raw.tsv
        gzip ${prefix}.raw.tsv
    fi

    echo "Done merging all results."
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$(echo \$(gzip --version 2>&1) | sed 's/^.*gzip //; s/ .*\$//')
    END_VERSIONS
    """
}
