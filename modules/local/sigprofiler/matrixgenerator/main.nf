process SIGPROFILER_MATRIXGENERATOR {
    tag "${task.ext.prefix}"
    label 'process_single'

    container 'docker.io/ferriolcalvet/sigprofiler:latest'

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
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def prefix = task.ext.prefix ?: "samples"
    """
    sigprofiler_matrix_generator.py \\
                ${prefix} \\
                ${params.vep_genome}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sigprofiler: 1.2.1
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        sigprofiler: 1.2.1
    END_VERSIONS
    """
}