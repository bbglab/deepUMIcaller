process SIGPROFILER_MATRIXGENERATOR {
    tag "${task.ext.prefix}"
    label 'process_single'

    container 'docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5'

    input:
    path (vcf)

    output:
    path("input_mutations/output/plots/*"), optional : true, emit: output_plots
    path("input_mutations/output/ID/*")   , optional : true, emit: matrices_ID
    path("input_mutations/output/DBS/*")  , optional : true, emit: matrices_DBS
    path("input_mutations/output/SBS/*")  , optional : true, emit: matrices_SBS
    path("input_mutations/output/TSB/*")  , optional : true, emit: transcription_bias
    path "versions.yml"                                    , topic: versions


    script:
    def prefix = task.ext.prefix ?: "samples"
    def args = task.ext.args ?: ""
    def genome = task.ext.genome_assembly ?: "GRCh38"
    """
    mkdir input_mutations
    cp *.vcf input_mutations/.

    SigProfilerMatrixGenerator matrix_generator \\
                ${prefix} \\
                ${genome} \\
                input_mutations/ \\
                ${args}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerMatrixGenerator: 1.3.5
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        SigProfilerMatrixGenerator: 1.3.5
    END_VERSIONS
    """
}