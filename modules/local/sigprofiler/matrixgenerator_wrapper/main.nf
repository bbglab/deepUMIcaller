// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process SIGPROFILER_MATRIXGENERATOR {
    tag "rrr"
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    // conda "YOUR-TOOL-HERE"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"
    container 'docker.io/ferriolcalvet/sigprofiler:human'

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    // tuple val(meta), path(vcf)
    path (vcf)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path("input_mutations/output/plots/*"), optional : true, emit: output_plots
    path("input_mutations/output/ID/*")   , optional : true, emit: matrices_ID
    path("input_mutations/output/DBS/*")  , optional : true, emit: matrices_DBS
    path("input_mutations/output/SBS/*")                   , emit: matrices_SBS
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix = task.ext.prefix ?: "samples"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    #!/usr/bin/env python3

    import os, sys, glob, gzip, shutil
    from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

    file_list = glob.glob(f"./*.vcf*")
    input_dir = "./input_mutations"

    # create dir
    os.mkdir(input_dir)
    
    print(file_list)

    for f in file_list:
        file = f.split("/")[-1]
        if file.endswith(".gz"):
            with gzip.open(f, 'r') as f_in, open(f'{input_dir}/{file[:-3]}', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        else:
            shutil.move(f, f'{input_dir}/{file}')

    matGen.SigProfilerMatrixGeneratorFunc("$prefix",
                                            "GRCh38",
                                            f"{input_dir}",
                                            plot=True,
                                            exome=False,
                                            bed_file=None,
                                            chrom_based=False,
                                            tsb_stat=False,
                                            seqInfo=False,
                                            cushion=100,
                                            $args
                                            )
    version_file = open("versions.yml", "w")
    version_file.write("${task.process}")
    version_file.write(f"    python: {sys.version.split(' ')[0]}")
    version_file.write("    sigprofiler: 1.2.1")
    version_file.close()
    """
}
