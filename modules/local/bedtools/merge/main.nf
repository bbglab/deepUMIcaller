process BEDTOOLS_MERGE {
    tag "$meta.id"
    label 'cpu_single'
    label 'time_low'
    label 'process_low_memory'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(vcf), path(bed)

    output:
    tuple val(meta), path('*.vcf_derived.bed')              , emit: vcf_bed
    tuple val(meta), path('*.vcf_derived.many.withID.bed')  , emit: vcf_bed_mut_ids
    tuple val(meta), path('*.regions_n_mutations.bed')      , emit: regions_plus_variants_bed
    path  "versions.yml"                                    , emit: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def amplify = task.ext.amplify ?: 5             // how many bases do you want to extend the region surrounding the variable

    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    //  awk '{ sum = length(\$4) + length(\$5); print \$1"\\t"\$2-$amplify"\\t"\$2 + sum }' | \\
    """
    grep -v '#' ${vcf} | \\
        awk '{ sum = length(\$4); print \$1"\\t"\$2-$amplify"\\t"\$2 + sum + $amplify"\\t"\$1";"\$2";"\$4";"\$5}' | \\
        sort -k1,1 -k2,3n \\
        > ${prefix}.vcf_derived.many.withID.bed
    
    cut -f-3 ${prefix}.vcf_derived.many.withID.bed > ${prefix}.vcf_derived.many.bed;

    bedtools \\
        merge \\
        -i ${prefix}.vcf_derived.many.bed \\
        -d 10 \\
        > ${prefix}.vcf_derived.bed;

    rm ${prefix}.vcf_derived.many.bed;

    bedtools \\
        merge \\
        -i <(cat ${bed} ${prefix}.vcf_derived.bed | cut -f -3 | sort -k1,1 -k2,3n) \\
        -d 10 \\
        ${args} \\
        > ${prefix}.regions_n_mutations.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}