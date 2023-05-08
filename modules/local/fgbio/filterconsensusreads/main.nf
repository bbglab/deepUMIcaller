process FGBIO_FILTERCONSENSUSREADS {
    tag "$meta.id"
    label 'process_medium_mem'

    conda "bioconda::fgbio=2.1.0 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fgbio:2.1.0--hdfd78af_0' :
        'biocontainers/fgbio:2.1.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(grouped_bam)
    path fasta

    output:
    tuple val(meta), path("*.filtered.bam")  , emit: bam
    path "versions.yml"                      , emit: versions

    script:
    def fgbio_args = task.ext.fgbio_args ?: ''
    def samtools_args = task.ext.samtools_args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mem_gb = 8
    if (!task.memory) {
        log.info '[fgbio FilterConsensusReads] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        mem_gb = task.memory.giga
    }
    // sort = false
    // if (sort) {
    //     fgbio_zipper_bams_output = "/dev/stdout"
    //     fgbio_zipper_bams_compression = 0 // do not compress if samtools is consuming it
    //     extra_command = " | samtools sort "
    //     extra_command += samtools_sort_args
    //     extra_command += " --template-coordinate"
    //     extra_command += " --threads "+ task.cpus
    //     extra_command += " -o " + prefix + ".filtered.bam##idx##"+ prefix + ".filtered.bam.bai"
    //     extra_command += " --write-index"
    //     extra_command += samtools_args
    // } else {
    fgbio_zipper_bams_output = prefix + ".filtered.bam"
    fgbio_zipper_bams_compression = 1
    extra_command = ""
    // }
    """
    fgbio \\
        -Xmx${mem_gb}g \\
        --tmp-dir=. \\
        --compression=${fgbio_zipper_bams_compression} \\
        FilterConsensusReads \\
        --input $grouped_bam \\
        --ref ${fasta} \\
        --output ${fgbio_zipper_bams_output} \\
        $fgbio_args;

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fgbio: \$( echo \$(fgbio --version 2>&1 | tr -d '[:cntrl:]' ) | sed -e 's/^.*Version: //;s/\\[.*\$//')
    END_VERSIONS
    """
}

// --max-base-error-rate ${max_base_error_rate} \\

// min-reads	M	Int	The minimum number of reads supporting a consensus base/read.	Required	3	 
// min-base-quality	N	PhredScore	Mask (make N) consensus bases with quality less than this threshold.	Required	1	 
// max-base-error-rate	e	Double	The maximum error rate for a single consensus base.	Required	3	0.1

// min-reads 3 1 1 
// min-base-quality 10
// max-base-error-rate 0.1

// min-reads	M	Int	The minimum number of reads supporting a consensus base/read.	Required	3	 
// min-base-quality	N	PhredScore	Mask (make N) consensus bases with quality less than this threshold.	Required	1	 
// max-read-error-rate	E	Double	The maximum raw-read error rate across the entire consensus read.	Required	3	0.025
// max-base-error-rate	e	Double	The maximum error rate for a single consensus base.	Required	3	0.1
// max-no-call-fraction	n	Double	Maximum fraction of no-calls in the read after filtering.	Optional	1	0.2
// min-mean-base-quality	q	PhredScore	The minimum mean base quality across the consensus read.	Optional	1	 
// require-single-strand-agreement	s	Boolean	Mask (make N) consensus bases where the AB and BA consensus reads disagree (for duplex-sequencing only).	Optional	1	false

// min-reads                3 1 1
// min-base-quality         10
// max-read-error-rate      0.025
// max-base-error-rate      0.1 0.3 0.3
// max-no-call-fraction     0.2
// min-mean-base-quality    20 // 
// require-single-strand-agreement	true


// --min-reads ${min_reads} \\
// --min-base-quality ${min_baseq} \\
// --max-base-error-rate 0.2 0.05 0.4 \\