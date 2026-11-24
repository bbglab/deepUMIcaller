process SPLITBAMCHROM {
    tag "$meta.id"
    label 'genomic_prep'
    
    conda "bioconda::samtools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
        'biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.chr*.bam"), emit: chrom_bams
    tuple val(meta), path("*.chr*.bam.bai"), emit: chrom_bais
    tuple val(meta), path("*.unknown.bam"), emit: unknown_bam
    tuple val(meta), path("*.unknown.bam.bai"), emit: unknown_bai
    tuple val(meta), path("*_filtering_stats.txt"), emit: filtering_stats
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Define chromosomes based on species
    def chromosomes = ""
    if (params.vep_species == "homo_sapiens") {
        chromosomes = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M"
    } else if (params.vep_species == "mus_musculus") {
        chromosomes = "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y M"
    } else {
        error "Unsupported species: ${params.vep_species}. Only 'homo_sapiens' and 'mus_musculus' are supported."
    }


    """
    # One pass over the index: contig, length, mapped, unmapped-placed.
    # Reused below for contig lists and read counts (no BAM decompression needed).
    samtools idxstats $bam > idxstats.txt
    cut -f1 idxstats.txt > contigs.txt

    # Total reads = mapped + unmapped (placed or not), as 'samtools view -c' would give
    total_reads=\$(awk '{s += \$3 + \$4} END {print s+0}' idxstats.txt)

    # Paired-but-not-proper-pair reads per contig, in a single pass over the BAM
    samtools view -@ $task.cpus -f 1 -F 2 $bam | cut -f3 | sort | uniq -c > unpaired_by_contig.txt

    # Helpers: read counts for the contigs listed in a file (\$1), from the input idxstats
    # or from the unpaired table; and read count of an indexed BAM (index only)
    count_idx()      { awk 'NR==FNR {want[\$1]=1; next} (\$1 in want) {s += \$3 + \$4} END {print s+0}' "\$1" idxstats.txt; }
    count_unpaired() { awk 'NR==FNR {want[\$1]=1; next} (\$2 in want) {s += \$1} END {print s+0}' "\$1" unpaired_by_contig.txt; }
    count_bam()      { samtools idxstats "\$1" | awk '{s += \$3 + \$4} END {print s+0}'; }

    echo "Chromosome Total_Reads Unpaired_Reads Paired_Reads Percentage_Unpaired((unpaired/unpaired+paired)*100) Percentage_Paired(paired/total)" > ${prefix}_filtering_stats.txt

    # Define chromosomes to process
    CHROMOSOMES='${chromosomes}'

    # Split for standard chromosomes and their associated sequences.
    # Outputs are intermediate: light compression (level=1) and index written in the same pass.
    for chrom in \$CHROMOSOMES; do
        grep -E '^chr'\$chrom'\$|^chr'\$chrom'_' contigs.txt > refs.txt || true
        if [ -s refs.txt ]; then
            samtools view $args -@ $task.cpus -b --output-fmt-option level=1 --write-index $bam \$(cat refs.txt) -o ${prefix}.chr\${chrom}.bam##idx##${prefix}.chr\${chrom}.bam.bai

            # Count reads for this chromosome
            chrom_total=\$(count_idx refs.txt)
            unpaired=\$(count_unpaired refs.txt)
            paired=\$(count_bam ${prefix}.chr\${chrom}.bam)
            percentage_unpaired=\$(awk -v unpaired="\$unpaired" -v paired="\$paired" 'BEGIN {total = unpaired + paired; if (total > 0) printf "%.2f", (unpaired / total) * 100; else print "0.00"}')
            percentage_paired=\$(awk -v paired="\$paired" -v total="\$chrom_total" 'BEGIN {if (total > 0) printf "%.2f", (paired / total) * 100; else print "0.00"}')
            echo "chr\${chrom} \$chrom_total \$unpaired \$paired \$percentage_unpaired \$percentage_paired" >> ${prefix}_filtering_stats.txt
        else
            echo "chr\${chrom} 0 0 0 0.00 0.00" >> ${prefix}_filtering_stats.txt
        fi
    done


    # Create unknown.bam with all other reads
    grep -E 'chrUn|HLA|chrEBV' contigs.txt > unknown_refs.txt || true
    if [ -s unknown_refs.txt ]; then
        # Create unknown.bam with chrUn, HLA, and chrEBV sequences (excluding unpaired)
        samtools view $args -@ $task.cpus -b --output-fmt-option level=1 --write-index $bam \$(cat unknown_refs.txt) -o ${prefix}.unknown.bam##idx##${prefix}.unknown.bam.bai

        # Count reads for unknown sequences
        unknown_total=\$(count_idx unknown_refs.txt)
        unpaired_unknown=\$(count_unpaired unknown_refs.txt)
        paired_unknown=\$(count_bam ${prefix}.unknown.bam)
        percentage_unpaired_unknown=\$(awk -v unpaired="\$unpaired_unknown" -v paired="\$paired_unknown" 'BEGIN {total = unpaired + paired; if (total > 0) printf "%.2f", (unpaired / total) * 100; else print "0.00"}')
        percentage_paired_unknown=\$(awk -v paired="\$paired_unknown" -v total="\$unknown_total" 'BEGIN {if (total > 0) printf "%.2f", (paired / total) * 100; else print "0.00"}')

        echo "unknown \$unknown_total \$unpaired_unknown \$paired_unknown \$percentage_unpaired_unknown \$percentage_paired_unknown" >> ${prefix}_filtering_stats.txt
    else
        touch ${prefix}.unknown.bam
        touch ${prefix}.unknown.bam.bai
        echo "unknown 0 0 0 0.00 0.00" >> ${prefix}_filtering_stats.txt
    fi

    #  Add total statistics
    total_unpaired=\$(awk '{sum += \$3} END {print sum}' ${prefix}_filtering_stats.txt)
    total_paired=\$(awk '{sum += \$4} END {print sum}' ${prefix}_filtering_stats.txt)
    percentage_total_unpaired=\$(awk -v unpaired="\$total_unpaired" -v paired="\$total_paired" 'BEGIN {total = unpaired + paired; if (total > 0) printf "%.2f", (unpaired / total) * 100; else print "0.00"}')
    percentage_total_paired=\$(awk -v paired="\$total_paired" -v total="\$total_reads" 'BEGIN {if (total > 0) printf "%.2f", (paired / total) * 100; else print "0.00"}')

    echo "Total \$total_reads \$total_unpaired \$total_paired \$percentage_total_unpaired \$percentage_total_paired" >> ${prefix}_filtering_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.chr1.bam ${prefix}.chr1.bam.bai ${prefix}.chr2.bam ${prefix}.chr2.bam.bai ${prefix}.chrX.bam ${prefix}.chrX.bam.bai ${prefix}.chrY.bam ${prefix}.chrY.bam.bai ${prefix}.unknown.bam ${prefix}.unknown.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}