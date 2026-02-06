process SPLIT_BED {
    tag "$meta.id"
    label 'bed_operations'

    input:
    tuple val(meta), path(targets_file)

    output:
    tuple val(meta), path("*_chunk_*"), emit: chunks
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "${meta.id}"
    def num_chunks = params.vardict_chunks
    """
    # Split the targets file into ${num_chunks} chunks
    total_lines=\$(wc -l < ${targets_file})
    
    # Calculate the number of lines per chunk
    lines_per_chunk=\$(( (total_lines + ${num_chunks} - 1) / ${num_chunks} ))
    
    # Split the file into chunks with the calculated number of lines
    # Using sample-specific prefix to avoid confusion in work directories
    split -l \${lines_per_chunk} ${targets_file} ${prefix}_chunk_
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coreutils: \$(echo \$(split --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
    END_VERSIONS
    """
}
