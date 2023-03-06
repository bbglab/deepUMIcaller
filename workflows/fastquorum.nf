/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFgcons.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.ref_fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }
if (params.ref_fasta) { ch_ref_fasta = Channel.fromPath(params.ref_fasta).collect() } else {
  log.error "No reference FASTA was specified (--ref_fasta)."
  exit 1
}

// The index directory is the directory that contains the FASTA
ch_ref_index_dir = ch_ref_fasta.map { it -> it.parent }

// Set various consensus calling and filtering parameters if not given
if (params.duplex_seq) {
  if (!params.groupreadsbyumi_strategy) { groupreadsbyumi_strategy = 'Paired' }
  else if (params.groupreadsbyumi_strategy != 'Paired') {
    log.error "config groupreadsbyumi_strategy must be 'Paired' for duplex-sequencing data"
    exit 1
  }
  if (!params.filter_min_reads) { filter_min_reads = '3 1 1' } else { filter_min_reads = params.filter_min_reads }
} else {
  if (!params.groupreadsbyumi_strategy) { groupreadsbyumi_strategy = 'Adjacency' }
  else if (params.groupreadsbyumi_strategy == 'Paired') {
    log.error "config groupreadsbyumi_strategy cannot be 'Paired' for non-duplex-sequencing data"
    exit 1
  } else {
	  groupreadsbyumi_strategy = params.groupreadsbyumi_strategy
  }
//   // here we should verify that both min reads parameters only have a single value inside a string
//   //     and are not for duplex reads
//   if (!params.call_min_reads IS A SINGLE NUMBER IN STRING FORMAT) {
//     log.error "config call_min_reads must be a single value in string format for non-duplex sequencing data"
//     exit 1 
//     }
//   if (!params.filter_min_reads IS A SINGLE NUMBER IN STRING FORMAT) {
//     log.error "config filter_min_reads must be a single value in string format for non-duplex sequencing data"
//     exit 1 
//     }
  if (!params.filter_min_reads) { filter_min_reads = '3' } else { filter_min_reads = params.filter_min_reads }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

include { ALIGN_BAM                         as ALIGNRAWBAM                 } from '../modules/local/align_bam_mod/main'
include { ALIGN_BAM                         as ALIGNCONSENSUSBAM           } from '../modules/local/align_bam_mod/main'
include { ALIGN_BAM                         as ALIGNDUPLEXCONSENSUSBAM     } from '../modules/local/align_bam_mod/main'

include { FGBIO_CLIPBAM                     as CLIPBAM                     } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMDUPLEX               } from '../modules/local/clipbam/main'

include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/local/fgbio/fastqtobam/main'

// include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI             } from '../modules/local/fgbio/groupreadsbyumi/main'
// include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/local/fgbio/groupreadsbyumi/main'

include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/local/fgbio/collectduplexseqmetrics/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSDUPLEX  } from '../modules/local/fgbio/filterconsensusreads/main'

include { CALLING_VARDICT                   as CALLINGVARDICT              } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTDUPLEX        } from '../modules/local/calling_vardict/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

include { FGBIO_SORTBAM                        as SORTBAM                  } from '../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_SORTBAM                        as SORTBAMCONS              } from '../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_SORTBAM                        as SORTBAMDUPLEXCONS        } from '../modules/nf-core/fgbio/sortbam/main'
include { FGBIO_SORTBAM                        as SORTBAMDUPLEXCONSFILT    } from '../modules/nf-core/fgbio/sortbam/main'
// include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/nf-core/fgbio/fastqtobam/main'

include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI             } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUSREADS    } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
// include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/nf-core/fgbio/filterconsensusreads/main'
// include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow FASTQUORUM {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )


    //
    // MODULE: Run fgbio FastqToBam
    //
    FASTQTOBAM(INPUT_CHECK.out.reads)
    // This is the unmapped BAM file: FASTQTOBAM.out.bam
    

    //
    // MODULE: Align with bwa mem
    //
    ALIGNRAWBAM(FASTQTOBAM.out.bam, ch_ref_index_dir, false)

    SORTBAM(ALIGNRAWBAM.out.bam)

    
    if (params.duplex_seq) {
        //
        // Run fgbio Duplex consensus pipeline
        //

        // MODULE: Run fgbio GroupReadsByUmi
        GROUPREADSBYUMIDUPLEX(SORTBAM.out.bam, "Paired")
        
        // MODULE: Run fgbio CallDuplexConsensusReads
        // CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam, call_min_reads, params.call_min_baseq)
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)

        // MODULE: Run fgbio CollecDuplexSeqMetrics
        COLLECTDUPLEXSEQMETRICS(GROUPREADSBYUMIDUPLEX.out.bam)
            
        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(CALLDUPLEXCONSENSUSREADS.out.bam, ch_ref_index_dir, false)   

        // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
        CLIPBAMDUPLEX(ALIGNDUPLEXCONSENSUSBAM.out.bam, ch_ref_fasta)

        //
        // ONLY DUPLEX READS
        //
        // MODULE: Run fgbio FilterConsensusReads
        FILTERCONSENSUSREADSDUPLEX(CLIPBAMDUPLEX.out.bam, ch_ref_fasta,
                                    filter_min_reads, params.filter_min_baseq,
                                    params.filter_max_base_error_rate)

        // MODULE: Sort BAM file
        SORTBAMDUPLEXCONSFILT(FILTERCONSENSUSREADSDUPLEX.out.bam)

        // Mutation calling for duplex reads
        CALLINGVARDICTDUPLEX(SORTBAMDUPLEXCONSFILT.out.bam, SORTBAMDUPLEXCONSFILT.out.index,
                            params.targetsfile,
                            ch_ref_fasta, ch_ref_index_dir)


        //
        // ALL READS
        //
        // MODULE: Sort BAM file
        SORTBAMDUPLEXCONS(CLIPBAMDUPLEX.out.bam)

        // Mutation calling for all reads
        CALLINGVARDICT(SORTBAMDUPLEXCONS.out.bam, SORTBAMDUPLEXCONS.out.index,
                        params.targetsfile,
                        ch_ref_fasta, ch_ref_index_dir)
        

    } else if (params.umi_only) {
        //
        // Run fgbio UMI-aware pipeline
        //

        // MODULE: Run fgbio GroupReadsByUmi
        GROUPREADSBYUMI(SORTBAM.out.bam, "Adjacency")

        // MODULE: Run fgbio CallMolecularConsensusReads
        // CALLMOLECULARCONSENSUSREADS(GROUPREADSBYUMI.out.bam, '1', params.call_min_baseq)
        CALLMOLECULARCONSENSUSREADS(GROUPREADSBYUMI.out.bam)

        // MODULE: Align with bwa mem
        ALIGNCONSENSUSBAM(CALLMOLECULARCONSENSUSREADS.out.bam, ch_ref_index_dir, false)

        // MODULE: Clip BAM file
        CLIPBAM(ALIGNCONSENSUSBAM.out.bam, ch_ref_fasta)
        
        // MODULE: Sort BAM file
        SORTBAMCONS(CLIPBAM.out.bam)

        // Mutation calling for non-duplex reads
        CALLINGVARDICT(SORTBAMCONS.out.bam, SORTBAMCONS.out.index,
                        params.targetsfile,
                        ch_ref_fasta, ch_ref_index_dir)
    }


    

    // Mutation calling for non-duplex reads
    // https://github.com/nf-core/sarek/blob/master/subworkflows/local/bam_variant_calling_tumor_only_mutect2/main.nf
    //  we could use this + force calling the mutations found in the duplex reads
    //      this needs to wait for the calling of duplex mutations to finish so that it can then know the positions to force call
    //      this should be easy to do just by passing as an argument the VCF file obtained from the duplex calling step.



    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFgcons.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(GROUPREADSBYUMIDUPLEX.out.histogram.map{it[1]}.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
