process CALLING_VARDICT_CHUNK {
    tag "$meta.id-${chunk_file.simpleName}"
    label 'variant_calling'
    
    conda "bioconda::vardict-java=1.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"    

    input:
    tuple val(meta), path(chunk_file), path(bam), path(bam_index)
    path fasta
    path fasta_dir

    output:
    tuple val(meta), path("*.genome.vcf"), emit: vcf
    tuple val(meta), path("*.raw.tsv")  , emit: raw
    path  "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def filter_args = task.ext.filter_args ?: ''
    def prefix = task.ext.prefix ?: ""
    // Use original sample ID without chunk suffix for prefix
    def sample_id = meta.id.replaceAll(/_chunk_[^_]+$/, '')
    prefix = "${sample_id}${prefix}"
    def chunk_name = chunk_file.simpleName
    """
    echo "Running vardict-java on chunk ${chunk_name}..."

    # Run vardict-java on this chunk
    vardict-java -G ${fasta_dir}/${fasta} \
        -N ${prefix} -b ${bam} \
        -c 1 -S 2 -E 3 -g 4 \
        $args \
        -th ${task.cpus} \
        ${chunk_file} > ${chunk_name}.raw.tsv

    echo "Vardict finished for chunk ${chunk_name}. Running teststrandbias..."

    # Process with teststrandbias or awk
    if ${params.use_teststrandbias}; then
        cat ${chunk_name}.raw.tsv | teststrandbias.R | var2vcf_valid.pl \
            -N ${prefix} $filter_args \
            > ${chunk_name}.genome.vcf
    else
        awk 'BEGIN { FS=OFS="\\t" } { 
            if (NF < 34) { 
                print "ERROR: Unexpected column count " NF > "/dev/stderr"; 
                exit 1; 
            } 
            for (i = 1; i <= 20; i++) printf "%s\\t", \$i; 
            printf "1.0\\t1.0\\t"; 
            for (i = 21; i <= NF; i++) { 
                printf "%s", \$i; 
                if (i < NF) printf "\\t";
            } 
            printf "\\n"; 
        }' ${chunk_name}.raw.tsv | var2vcf_valid.pl \
            -N ${prefix} -A -E -f 0.0 -p 0 -m 20 -v 2 \
            > ${chunk_name}.genome.vcf
    fi

    echo "Done processing chunk ${chunk_name}."
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$(vardict-java -h 2>&1 | grep -oP 'VarDict version \\K[0-9.]+' || echo "1.8.3")
    END_VERSIONS
    """
}
