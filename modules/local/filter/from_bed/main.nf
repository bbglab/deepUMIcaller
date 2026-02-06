process FILTER_FROM_BED {
    tag "$meta.id"
    label 'variant_filtering'

    // Use either conda or container, depending on profile
    conda "bioconda::pybedtools=0.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' : 
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"

    input:
    tuple val(meta), path(vcf_file), path(vcf_derived_bed), path(mask_bed), val(filter_name)

    output:
    tuple val(meta), path("${meta.id}${task.ext.prefix ?: ''}.${filter_name}.vcf"), path(vcf_derived_bed)    , emit: filtered_vcf_bed
    tuple val(meta), path("${meta.id}${task.ext.prefix ?: ''}.${filter_name}.vcf")                           , emit: filtered_vcf
    path  "versions.yml"                                                                                     , topic: versions

    script:
    def prefix = task.ext.prefix ?: ""
    def outfile = "${meta.id}${prefix}.${filter_name}.vcf"
    def bedfile = "${meta.id}${prefix}.${filter_name}_file.bed"
    """
    if [[ "${mask_bed}" == *.gz ]]; then
        bedtools intersect -a ${vcf_derived_bed} -b <(zcat ${mask_bed}) -u > ${bedfile}
    else
        bedtools intersect -a ${vcf_derived_bed} -b ${mask_bed} -u > ${bedfile}
    fi

    # If there is nothing in the intersection, just copy the VCF
    if [ -s ${bedfile} ]; then
        add_filter_from_bed.py \\
            ${vcf_file} \\
            ${bedfile} \\
            ${outfile} \\
            ${filter_name}
    else
        cp ${vcf_file} ${outfile}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    def outfile = "${meta.id}${prefix}.${filter_name}.vcf"
    """
    touch ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
