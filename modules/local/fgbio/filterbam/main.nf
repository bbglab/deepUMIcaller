process FGBIO_FILTERBAM {
    tag "$meta.id"
    label 'process_medium_mem'

    // TODO
    // update with only fgbio and samtools
    conda "bioconda::fgbio=2.1.0 bioconda::bwa=0.7.17 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ? 
        'https://depot.galaxyproject.org/singularity/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' : 
        'biocontainers/mulled-v2-69f5207f538e4de9ef3bae6f9a95c5af56a88ab8:82d3ec41f9f1227f7183d344be46f73365efa704-0' }"

    input:
    tuple val(meta), path(bam), path (index)
    path regions

    output:
    tuple val(meta), path("*.txt"), emit: read_names
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def samtools_args = task.ext.samtools_args ?: ''
    def fgbio_args = task.ext.fgbio_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fgbio_mem_gb = 4

    if (!task.memory) {
        log.info '[fgbio FilterBam] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else if (fgbio_mem_gb > task.memory.giga) {
        if (task.memory.giga < 2) {
            fgbio_mem_gb = 1
        } else {
            fgbio_mem_gb = task.memory.giga - 1
        }
    }
    fgbio_compression = 0
    """
    fgbio -Xmx${fgbio_mem_gb}g \\
        --compression ${fgbio_compression} \\
        --async-io=true \\
        FilterBam \\
        --input ${bam} \\
        --intervals ${regions} \\
        --output /dev/stdout \\
        --remove-duplicates false \\
        ${fgbio_args} \\
        | samtools view ${samtools_args} - | cut -f1 > ${prefix}.read_ids.txt
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
// TODO
// add final step to the process in order to remove the duplicated read names
// rev file.txt | sort -u | rev > file_unique.txt
// the rev trick here allows a faster sorting and duplicate removal by starting to compare the reads from the end
//  this is more efficient since the beginning of the read name is shared across many more reads than the end

// intervals	l	PathToIntervals	Optionally remove reads not overlapping intervals.	Optional	1
// remove-duplicates	D	Boolean	If true remove all reads that are marked as duplicates.	Optional	1	true

// remove-unmapped-reads	U	Boolean	Remove all unmapped reads.	Optional	1	true
// min-map-q	M	Int	Remove all mapped reads with MAPQ lower than this number.	Optional	1	1
// remove-single-end-mappings	P	Boolean	Removes non-PE reads and any read whose mate pair is unmapped.	Optional	1	false
// remove-secondary-alignments	S	Boolean	Remove all reads marked as secondary alignments.	Optional	1	true
// min-insert-size	 	Int	Remove all reads with insert size < this value.	Optional	1	 
// max-insert-size	 	Int	Remove all reads with insert size > this value.	Optional	1	 
// min-mapped-bases	m	Int	Remove reads with fewer than this many mapped bases.	Optional	1	 