process NS_X_POSITION {
    tag "$meta.id"
    label 'cpu_single'
    label 'time_low'
    label 'memory_medium'

    conda "bioconda::tabix=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tabix:1.11--hdfd78af_0' :
        'biocontainers/tabix:1.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(pileup), path(pileuptabix)

    output:
    tuple val(meta), path("*.tsv.gz"), path("*.tsv.gz.tbi") , emit: ns_tsv
    path  "versions.yml"                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    zcat $pileup | \\
            awk 'BEGIN{FS=OFS="\\t"} {print \$1"\\t"\$2"\\t"\$4"\\t"gsub(/N/,"",\$5)"\\t"gsub(/n/,"",\$5)"\\t"gsub(/\\*/,"",\$5)}' | \\
            awk '{\$3=\$3-\$6;\$4=\$4+\$5; print \$1"\\t"\$2"\\t"\$3"\\t"\$4}' | \\
            bgzip \\
            > ${prefix}.Ns_per_position.tsv.gz;
    tabix -s 1 -b 2 -e 2 ${prefix}.Ns_per_position.tsv.gz;
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Ns_per_position.tsv.gz
    touch ${prefix}.Ns_per_position.tsv.gz.tbi;

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
