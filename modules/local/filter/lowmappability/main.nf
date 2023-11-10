// /workspace/projects/bladder_ts/scripts/add_filters_vcf

process FILTER_LOW_MAPPABILITY {
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

    output:
    tuple val(meta), path("*.low_mappable.vcf"), emit: filtered_vcf
    path  "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def low_mappable_bed = task.ext.low_mappable ?: ''

    // TODO add condition to check if the file argument is empty, if it is, return same vcf file
    // TODO change it so that we can provide a BED file of the unmappable genome, not having to provide the mappable one
    """
    bedtools intersect -a ${vcf_derived_bed} -b ${low_mappable_bed} -v  > ${prefix}.lowmappable_file.bed

    # if there is nothing in the intersection do not filter the VCF file
    if [ -s ${prefix}.lowmappable_file.bed ]; then
        add_filter_lowmappability.py \\
                ${vcf_file} \\
                ${prefix}.lowmappable_file.bed \\
                ${prefix}.low_mappable.vcf \\
                low_mappability;
    else
        cp ${vcf_file} ${prefix}.low_mappable.vcf;
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.low_mappable.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}

