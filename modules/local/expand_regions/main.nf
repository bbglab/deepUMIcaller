process EXPAND_PANEL {
    tag "$meta.id"
    label 'postprocess_compute'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(bedfile)

    output:
    tuple val(meta), path("*expanded.bed") , emit: expanded_bed
    path  "versions.yml"                   , topic: versions


    script:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    awk -v OFS='\\t' '{\$2=\$2-250;\$3=\$3+250;print \$0}' ${bedfile} > ${prefix}.expanded250.bed6.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    touch ${prefix}.expanded250.bed6.bed;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}

// tabix -s 1 -b 2 -e 2 ${prefix}.Ns_per_position.tsv.gz;
// tabix: \$(echo \$(tabix --version 2>&1) | sed 's/^.*tabix (htslib) //; s/Copyright.*\$//')
// zcat $pileup | \\
//         awk 'BEGIN{FS=OFS="\\t"} {print \$1"\\t"\$2"\\t"gsub(/N/,"",\$5)"\\t"gsub(/n/,"",\$5)}' | \\
//         awk '{\$3=\$3+\$4; print \$1"\\t"\$2"\\t"\$3}' | \\
//         bgzip | \\
//         > ${prefix}.Ns_per_position.tsv.gz

// bedtools intersect -u -a <(zcat $MPILEUP_PATH/$SAMPLE.mpileup_out.tsv.gz |
//                             awk 'BEGIN{FS=OFS="\t"} {print $1"\t"$2"\t"$4"\t"gsub(/N/,"",$5)"\t"gsub(/n/,"",$5)}' |
//                             awk '{$4=$4+$5; print $1"\t"$2"\t"$2"\t"$3"\t"$4}') -b $BED_FILE |
//                             gzip > $N_DEPTH_PATH/$SAMPLE.high.Ns_per_position.tsv.gz
