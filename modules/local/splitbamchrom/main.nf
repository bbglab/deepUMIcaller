process SPLITBAMCHROM {
    tag "$meta.id"
    label 'process_medium'
    
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

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Get total reads in the input BAM
    total_reads=\$(samtools view -@ $task.cpus -c $bam)
    echo "Chromosome Total_Reads Unpaired_Reads Paired_Reads Percentage_Unpaired((unpaired/unpaired+paired)*100) Percentage_Paired(paired/total)" > ${prefix}_filtering_stats.txt

    # Define chromosomes to process
    CHROMOSOMES='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M'

    # Split for standard chromosomes and their associated sequences
    for chrom in \$CHROMOSOMES; do
        refs=\$(samtools idxstats $bam | cut -f1 | grep -E '^chr'\$chrom'\$|^chr'\$chrom'_')
        if [ ! -z "\$refs" ]; then
            samtools view  -f 3 -F 2304  $args -@ $task.cpus -b $bam \$refs -o ${prefix}.chr\${chrom}.bam
            samtools index -@ $task.cpus ${prefix}.chr\${chrom}.bam

            # Count reads for this chromosome
            chrom_total=\$(samtools view -@ $task.cpus -c $bam \$refs)
            unpaired=\$(samtools view -@ $task.cpus -c -f 1 -F 2 $bam \$refs)
            paired=\$(samtools view -@ $task.cpus -c -f 3 -F 2304 ${prefix}.chr\${chrom}.bam)
            percentage_unpaired=\$(awk -v unpaired="\$unpaired" -v paired="\$paired" 'BEGIN {total = unpaired + paired; if (total > 0) printf "%.2f", (unpaired / total) * 100; else print "0.00"}')
            percentage_paired=\$(awk -v paired="\$paired" -v total="\$chrom_total" 'BEGIN {if (total > 0) printf "%.2f", (paired / total) * 100; else print "0.00"}')
            echo "chr\${chrom} \$chrom_total \$unpaired \$paired \$percentage_unpaired \$percentage_paired" >> ${prefix}_filtering_stats.txt
        else
            echo "chr\${chrom} 0 0 0 0.00 0.00" >> ${prefix}_filtering_stats.txt
        fi
    done


    # Create unknown.bam with all other reads
    samtools idxstats $bam | cut -f1 | grep -E 'chrUn|HLA|chrEBV' > unknown_refs.txt
    if [ -s unknown_refs.txt ]; then
        # Create unknown.bam with chrUn, HLA, and chrEBV sequences (excluding unpaired)
        samtools view  -f 3 -F 2304  $args -@ $task.cpus -b $bam \$(cat unknown_refs.txt) -o ${prefix}.unknown.bam
        samtools index -@ $task.cpus ${prefix}.unknown.bam

        # Count reads for unknown sequences
        unknown_total=\$(samtools view -@ $task.cpus -c $bam \$(cat unknown_refs.txt))
        unpaired_unknown=\$(samtools view -@ $task.cpus -c -f 1 -F 2 $bam \$(cat unknown_refs.txt))
        paired_unknown=\$(samtools view -@ $task.cpus -c ${prefix}.unknown.bam)
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