// /workspace/projects/bladder_ts/scripts/add_filters_vcf

process FILTER_LOW_COMPLEXITY {
    tag "$meta.id"
    label 'cpu_single'
    label 'time_low'
    label 'process_low_memory'
    
    conda "bioconda::pybedtools=0.9.1--py38he0f268d_0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
            'https://depot.galaxyproject.org/singularity/pybedtools:0.9.1--py38he0f268d_0' : 
            'biocontainers/pybedtools:0.9.1--py38he0f268d_0' }"


    input:
    tuple val(meta), path(vcf_file), path(vcf_derived_bed)
    path (low_complex_bed)

    output:
    tuple val(meta), path("*.low_complex.vcf"), path(vcf_derived_bed)   , emit: filtered_vcf_bed
    path  "versions.yml"                                                , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    bedtools intersect -a ${vcf_derived_bed} -b ${low_complex_bed} -u  > ${prefix}.lowcomplexrep_file.bed

    # if there is nothing in the intersection do not filter the VCF file
    if [ -s ${prefix}.lowcomplexrep_file.bed ]; then
        add_filter_lowcomplexrep.py \\
                        ${vcf_file} \\
                        ${prefix}.lowcomplexrep_file.bed \\
                        ${prefix}.low_complex.vcf \\
                        low_complex_repetitive;
    else
        cp ${vcf_file} ${prefix}.low_complex.vcf;
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.low_complex.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

