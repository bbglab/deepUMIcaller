process FGUMI_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'consensus_filter'

    conda "bioconda::fgumi"
    label 'fgumi_tools' 

    input:
    tuple val(meta), path(grouped_bam)
    path fasta

    output:
    tuple val(meta), path("*.filtered.bam")  , emit: bam
    path "versions.yml"                      , topic: versions

    script:
    def fgumi_args = task.ext.fgumi_args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    def Ns_per_read = (params.maxN_per_read + params.left_clip + params.right_clip) / params.read_length_minus_tag
    def max_no_call_fraction = Ns_per_read ? "--max-no-call-fraction ${Ns_per_read}" : "--max-no-call-fraction 0.2"
    """
    fgumi filter \\
        --input $grouped_bam \\
        --ref ${fasta} \\
        --output ${prefix}.filtered.bam \\
        --threads ${task.cpus} \\
        ${max_no_call_fraction} \
        ${fgumi_args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgumi: \$(fgumi --version | sed 's/^fgumi //')
    END_VERSIONS
    """
}
