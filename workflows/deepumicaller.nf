/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowFgcons.initialise(params, log)

// TODO nf-core:
// Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.ref_fasta, params.targetsfile  ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Check mandatory parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.input) { ch_input = file(params.input) } else {
    exit 1, 'Input samplesheet not specified!'
    }

if (params.ref_fasta) {
    ch_ref_fasta = Channel.fromPath(params.ref_fasta).collect()

    // define additional fasta file names
    ch_ref_fasta_file = file(params.ref_fasta, checkIfExists: true)
    ch_ref_fasta_fai_index = file("${ch_ref_fasta_file}.fai", checkIfExists: true)
    ch_ref_fasta_dict = file("${ch_ref_fasta_file.parent/ch_ref_fasta_file.baseName}.dict", checkIfExists: true)

} else {
    log.error "No reference FASTA was specified (--ref_fasta)."
    exit 1
    }


// The index directory is the directory that contains the FASTA
ch_ref_index_dir = ch_ref_fasta.map { it -> it.parent }
// TODO
// check if the index file for the reference genome is present
// if (ch_ref_index_dir) { file("${file(params.ref_fasta).parent}/${file(params.ref_fasta).name}.amb", checkIfExists: true) }



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
include { INPUT_CHECK                                                      } from '../subworkflows/local/input_check'

include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/local/fgbio/fastqtobam/main'

include { ALIGN_BAM                         as ALIGNRAWBAM                 } from '../modules/local/align_bam/main'

include { ALIGN_BAM                         as ALIGNCONSENSUSBAM           } from '../modules/local/align_bam/main'
include { ALIGN_BAM                         as ALIGNDUPLEXCONSENSUSBAM     } from '../modules/local/align_bam/main'

include { FGBIO_FILTERBAM                   as FGSELECTREADS               } from '../modules/local/fgbio/filterbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAM                     } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMLOWCONF              } from '../modules/local/clipbam/main'

include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/local/fgbio/collectduplexseqmetrics/main'

// include { PLOTDUPLEXMETRICS                 as PLOTDUPLEXMETRICS           } from '../modules/local/duplexfamilymetrics/main'
// include { FAMILYMETRICS                     as FAMILYMETRICS               } from '../modules/local/familymetrics/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSDUPLEX  } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSLOWCONF } from '../modules/local/fgbio/filterconsensusreads/main'

include { CALLING_VARDICT                   as CALLINGVARDICT              } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTDUPLEX        } from '../modules/local/calling_vardict/main'

// include { BBGPOSTANALYSIS                     as BBGPOSTANALYSIS               } from '../modules/local/BBGPOSTANALYSIS/main'

include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTDUPLEX           } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOT                 } from '../modules/local/sigprofiler/matrixgenerator/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTP                             as FASTP                       } from '../modules/nf-core/fastp/main'

include { FASTQC                            as FASTQC                      } from '../modules/nf-core/fastqc/main'
// include { FASTQC                            as POSTTRIMQC                  } from '../modules/nf-core/fastqc/main'

include { BAMUTIL_TRIMBAM                   as TRIMBAM                     } from '../modules/nf-core/bamutil/trimbam/main'

include { PICARD_BEDTOINTERVALLIST          as BEDTOINTERVAL               } from '../modules/nf-core/picard/bedtointervallist/main'

//  Metrics
include { QUALIMAP_BAMQC                    as QUALIMAPQC                  } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCDUPLEX            } from '../modules/nf-core/qualimap/bamqc/main'
include { SAMTOOLS_DEPTH                    as COMPUTEDEPTH                } from '../modules/nf-core/samtools/depth/main'
include { PICARD_COLLECTMULTIPLEMETRICS     as COLLECTMULTIPLEMETRICS      } from '../modules/nf-core/picard/collectmultiplemetrics/main'

// Versions and reports
include { MULTIQC                                                          } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                      } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Sorting
include { SAMTOOLS_SORT                     as SORTBAM                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMFIXED                } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCONS                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONS           } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSFILT       } from '../modules/nf-core/samtools/sort/main'

include { SAMTOOLS_VIEW                     as FILTERBAM                   } from '../modules/nf-core/samtools/view/main'

// include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/nf-core/fgbio/fastqtobam/main'

include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI             } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUSREADS    } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
// include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/nf-core/fgbio/filterconsensusreads/main'
// include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'


// Download annotation cache if needed
include { PREPARE_CACHE                                                   } from '../subworkflows/local/prepare_cache/main'


// Annotation
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEALL             } from '../subworkflows/local/vcf_annotate_all/main'
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEALLDUPLEX       } from '../subworkflows/local/vcf_annotate_all/main'









/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DEEPUMICALLER {

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // Download Ensembl VEP cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = params.vep_cache ? [] : Channel.of([ [ id:"${params.vep_genome}.${params.vep_cache_version}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    if (params.download_cache) {
        PREPARE_CACHE(ensemblvep_info)
        vep_cache = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }
        ch_versions = ch_versions.mix(PREPARE_CACHE.out.versions)
    } else {
        vep_cache = params.vep_cache
    }
    vep_extra_files = []
    
    
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    // READ PREPROCESSING
    if (params.trim_adapters){
        // MODULE: Run FASTP
        FASTP(INPUT_CHECK.out.reads,
                        [], // we are not using any adapter fastas at the moment
                        false,
                        false)
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        reads_to_qc = FASTP.out.reads
    } else {
        reads_to_qc = INPUT_CHECK.out.reads
    }


    // MODULE: Run FastQC
    FASTQC (
        reads_to_qc
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))


    // MODULE: Run fgbio FastqToBam
    // to get the UMIs out of the reads and into the tag
    FASTQTOBAM(reads_to_qc)
    ch_versions = ch_versions.mix(FASTQTOBAM.out.versions.first())
    // This is the unmapped BAM file: FASTQTOBAM.out.bam


    // Decide whether we clip the beginning and/or end of the reads or nothing
    if ( (params.left_clip > 0) || (params.right_clip > 0) ) {
        // params.left_clip = 4     ->      Remove 4bp from the 5' end of the reads
        TRIMBAM(FASTQTOBAM.out.bam, params.left_clip, params.right_clip)
        ch_versions = ch_versions.mix(TRIMBAM.out.versions.first())
        bam_to_align = TRIMBAM.out.bam
    } else {
        bam_to_align = FASTQTOBAM.out.bam
    }


    // MODULE: Align with bwa mem
    ALIGNRAWBAM(bam_to_align, ch_ref_index_dir, false)
    ch_versions = ch_versions.mix(ALIGNRAWBAM.out.versions.first())

    SORTBAM(ALIGNRAWBAM.out.bam)
    ch_versions = ch_versions.mix(SORTBAM.out.versions.first())

    // COLLECTMULTIPLEMETRICS(SORTBAM.out.bam, SORTBAM.out.csi.map{it -> it [1]}, ch_ref_fasta, ch_ref_fasta_fai_index)
    // ch_versions = ch_versions.mix(COLLECTMULTIPLEMETRICS.out.versions.first())

    if (params.targetsfile){
        QUALIMAPQC(SORTBAM.out.bam, params.targetsfile)
        ch_versions = ch_versions.mix(QUALIMAPQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())

        targets_bed = Channel.of([ [ id:"${file(params.targetsfile).getSimpleName()}" ], file(params.targetsfile) ])
        BEDTOINTERVAL(targets_bed, ch_ref_fasta_dict, [])
        ch_versions = ch_versions.mix(BEDTOINTERVAL.out.versions.first())

        // TODO
        // collect metrics for on target vs off target reads

        // truncate BAM to keep only the reads that are on target
        FGSELECTREADS(SORTBAM.out.bam, BEDTOINTERVAL.out.interval_list.first().map{it -> it [1]})
        FILTERBAM(SORTBAM.out.bam, SORTBAM.out.csi.map{it -> it[1]},  [], FGSELECTREADS.out.read_names.map{it -> it[1]} )
        SORTBAMFIXED(FILTERBAM.out.bam)

        bam_to_group = SORTBAMFIXED.out.bam

    } else {
        QUALIMAPQC(SORTBAM.out.bam, [])
        ch_versions = ch_versions.mix(QUALIMAPQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())

        bam_to_group = SORTBAM.out.bam

    }



    if (params.duplex_seq) {
        //
        // Run fgbio Duplex consensus pipeline
        //

        // MODULE: Run fgbio GroupReadsByUmi
        GROUPREADSBYUMIDUPLEX(bam_to_group, "Paired")
        ch_versions = ch_versions.mix(GROUPREADSBYUMIDUPLEX.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(GROUPREADSBYUMIDUPLEX.out.histogram.map{it[1]}.collect())


        // MODULE: Run fgbio CallDuplexConsensusReads
        // CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam, call_min_reads, params.call_min_baseq)
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)
        ch_versions = ch_versions.mix(CALLDUPLEXCONSENSUSREADS.out.versions.first())

        // MODULE: Run fgbio CollecDuplexSeqMetrics
        COLLECTDUPLEXSEQMETRICS(GROUPREADSBYUMIDUPLEX.out.bam)
        ch_versions = ch_versions.mix(COLLECTDUPLEXSEQMETRICS.out.versions.first())

        // TODO
        // add metrics plots module
        // COLLECTDUPLEXSEQMETRICS.out.metrics // the problem here is that there are many files, we only need one
        // COLLECTDUPLEXSEQMETRICS.out.specific_metrics

        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(CALLDUPLEXCONSENSUSREADS.out.bam, ch_ref_index_dir, false)

        //
        // ONLY DUPLEX READS
        //
        // MODULE: Run fgbio FilterConsensusReads
        FILTERCONSENSUSREADSDUPLEX(ALIGNDUPLEXCONSENSUSBAM.out.bam, ch_ref_fasta)
        ch_versions = ch_versions.mix(FILTERCONSENSUSREADSDUPLEX.out.versions.first())

        // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
        CLIPBAM(FILTERCONSENSUSREADSDUPLEX.out.bam, ch_ref_fasta)
        ch_versions = ch_versions.mix(CLIPBAM.out.versions.first())

        // MODULE: Sort BAM file
        SORTBAMDUPLEXCONSFILT(CLIPBAM.out.bam)

        // Compute depth of the "raw" reads aligned to the genome
        COMPUTEDEPTH(SORTBAMDUPLEXCONSFILT.out.bam)

        // Quality check
        QUALIMAPQCDUPLEX(SORTBAMDUPLEXCONSFILT.out.bam, params.targetsfile)
        ch_versions = ch_versions.mix(QUALIMAPQCDUPLEX.out.versions.first())
        // ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCDUPLEX.out.results.map{it[1]}.collect())
    

        // Mutation calling for duplex reads
        CALLINGVARDICTDUPLEX(SORTBAMDUPLEXCONSFILT.out.bam, SORTBAMDUPLEXCONSFILT.out.csi,
                            params.targetsfile,
                            ch_ref_fasta, ch_ref_index_dir)
        ch_versions = ch_versions.mix(CALLINGVARDICTDUPLEX.out.versions.first())


        VCFANNOTATEALLDUPLEX(CALLINGVARDICTDUPLEX.out.vcf,
                            ch_ref_fasta,
                            "GRCh38",
                            "homo_sapiens", 
                            "108",
                            vep_cache,
                            vep_extra_files)
        ch_versions = ch_versions.mix(VCFANNOTATEALLDUPLEX.out.versions.first())

        CALLINGVARDICTDUPLEX.out.vcf.map{it -> it[1]}.set { mutation_files_duplex }

        SIGPROFPLOTDUPLEX(mutation_files_duplex.collect())
        // ch_versions = ch_versions.mix(SIGPROFPLOTDUPLEX.out.versions.first())

        // TODO
        // add a module that plots the signature, with the correct normalization
        // use bgsignature, and the input targets file


        //
        // ALL READS
        //
        if (params.duplex_low_conf){

            // TODO
            // choose filtering step parameters
            FILTERCONSENSUSREADSLOWCONF(ALIGNDUPLEXCONSENSUSBAM.out.bam, ch_ref_fasta)

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMLOWCONF(FILTERCONSENSUSREADSLOWCONF.out.bam, ch_ref_fasta)
            ch_versions = ch_versions.mix(CLIPBAMLOWCONF.out.versions.first())

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONS(CLIPBAMLOWCONF.out.bam)

            // TODO
            // decide whether we want to add a step to compute metics on the BAM file, either qualimap or depth

            // Mutation calling for all reads
            CALLINGVARDICT(SORTBAMDUPLEXCONS.out.bam, SORTBAMDUPLEXCONS.out.csi,
                            params.targetsfile,
                            ch_ref_fasta, ch_ref_index_dir)

            VCFANNOTATEALL(CALLINGVARDICT.out.vcf,
                            ch_ref_fasta,
                            "GRCh38",
                            "homo_sapiens", 
                            "108",
                            vep_cache,
                            vep_extra_files)

            CALLINGVARDICT.out.vcf.map{it -> it[1]}.set { mutation_files }
            SIGPROFPLOT(mutation_files.collect())
        }

        // TODO
        // FGBIO FILTER SOMATIC VARIANTS
        // http://fulcrumgenomics.github.io/fgbio/tools/latest/FilterSomaticVcf.html


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
        CALLINGVARDICT(SORTBAMCONS.out.bam, SORTBAMCONS.out.csi,
                        params.targetsfile,
                        ch_ref_fasta, ch_ref_index_dir)
    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowFgcons.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

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
