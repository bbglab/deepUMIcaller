//
// Subworkflow with functionality specific to the bbglab/intogenplus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )


    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text =  logo(monochrome_logs)
    post_help_text = '\n' + citation(workflow) + '\n'
    def workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../bbglab> --input path/to/input/ --outdir path/to/output/"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )


    emit:
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}



def logo(monochrome_logs) {
    String.format(
        """
        ${monochrome_logs ? '' : '\033[0;33m'}
        888                             888     888 888b     d888 888                   888 888
        888                             888     888 8888b   d8888 888                   888 888
        888                             888     888 88888b.d88888 888                   888 888
    .d88888  .d88b.   .d88b.  88888b.   888     888 888Y88888P888 888  .d8888b  8888b.  888 888  .d88b.  888d888
    d88" 888 d8P  Y8b d8P  Y8b 888 "88b  888     888 888 Y888P 888 888 d88P"        "88b 888 888 d8P  Y8b 888P"
    888  888 88888888 88888888 888  888  888     888 888  Y8P  888 888 888      .d888888 888 888 88888888 888
    Y88b 888 Y8b.     Y8b.     888 d88P  Y88b. .d88P 888   "   888 888 Y88b.    888  888 888 888 Y8b.     888
    "Y88888  "Y8888   "Y8888  88888P"    "Y88888P"  888       888 888  "Y8888P "Y888888 888 888  "Y8888  888
                            888
                            888
                            888
        ${monochrome_logs ? '' : '\033[0m'}
        """.stripIndent()
    )
}


def citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
}
