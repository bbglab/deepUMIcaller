process MUTS_PER_POS {
    tag "$meta.id"
    label 'time_medium'
    label 'cpu_low'
    label 'process_medium_high_memory'

    container 'docker.io/ferriolcalvet/pysam'

    input:
    tuple val(meta), path(bam), path(bam_index), path(vcf)

    output:
    tuple val(meta), path("**.png")                         , emit: all_plots
    tuple val(meta), path("**_BasePerPosWithoutNs_mqc.png") , emit: plots
    tuple val(meta), path("**_MutsPerCycle_mqc.csv")     , emit: positions_csv
    tuple val(meta), path("**")                          , emit: others
    path  "versions.yml"                                 , topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    grep -v '##' $vcf > ${prefix}.no_header.vcf;
    count_muts_per_cycle.py \\
                --inFile ${bam} \\
                --inVCF ${prefix}.no_header.vcf \\
                -o ${prefix} \\
                ${args}
    
    # Rename files for MultiQC recognition
    mv ${prefix}_MutsPerCycle.dat.csv ${prefix}_MutsPerCycle_mqc.csv
    mv ${prefix}_BasePerPosWithoutNs.png ${prefix}_BasePerPosWithoutNs_mqc.png
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}_BasePerPosInclNs_mqc.png
    touch ${prefix}_BasePerPosWithoutNs_mqc.png
    touch ${prefix}_MutsPerCycle_mqc.csv
    touch ${prefix}_mutsPerRead_mqc.png
    touch ${prefix}.no_header.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

// vcf is filtered for only low VAF variants in principle,
//      but we could let the python script filter it itself

// this was run here in the past
// /data/bbg/projects/prominent/analysis/dev_pipeline/prom10/get_muts_per_position
