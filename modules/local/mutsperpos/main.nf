process MUTS_PER_POS {
    tag "$meta.id"
    label 'time_medium'
    label 'cpu_low'
    label 'process_medium_high_memory'
    // conda "bioconda::pysam=0.21.0--py38h15b938a_1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/pysam:0.21.0--py38h15b938a_1' :
    //     'biocontainers/pysam:0.21.0--py38h15b938a_1' }"
    container 'docker.io/ferriolcalvet/pysam'

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf)

    output:
    tuple val(meta), path("**.png")     , emit: plots
    tuple val(meta), path("**")         , emit: others
    path  "versions.yml"                , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_size = task.ext.read_size ?: "142"
    """
    # vcf is filtered for only low VAF variants in principle, but we could let the python script filter it itself
    grep -v '##' $vcf > ${prefix}.no_header.vcf;
    count_muts_per_cycle.py \\
                --inFile ${bam} \\
                --inVCF ${prefix}.no_header.vcf \\
                -o ${prefix} \\
                -l ${read_size} \\
                --filter INCLUDE \\
                -t 0 \\
                --clonality_limit 0.1 \\
                -c
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_BasePerPosInclNs.png
    touch ${prefix}_BasePerPosWithoutNs.png
    touch ${prefix}_MutsPerCycle.dat.csv
    touch ${prefix}.mutsPerRead.png
    touch ${prefix}.no_header.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

// this was run here in the past
// /workspace/projects/prominent/analysis/dev_pipeline/prom10/get_muts_per_position
