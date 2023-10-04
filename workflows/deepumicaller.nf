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

include { BAM_FILTER_READS                                                 } from '../subworkflows/local/bam_filter_reads/main'


include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/local/fgbio/fastqtobam/main'

include { ALIGN_BAM                         as ALIGNRAWBAM                 } from '../modules/local/align_bam/main'

include { ALIGN_BAM                         as ALIGNCONSENSUSBAM           } from '../modules/local/align_bam/main'
include { ALIGN_BAM                         as ALIGNDUPLEXCONSENSUSBAM     } from '../modules/local/align_bam/main'

include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/local/fgbio/collectduplexseqmetrics/main'
include { FAMILYSIZEMETRICS                 as FAMILYMETRICS               } from '../modules/local/familymetrics/main'

include { FGBIO_CLIPBAM                     as CLIPBAMLOW                  } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMMED                  } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMHIGH                 } from '../modules/local/clipbam/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSLOW     } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSMED     } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSHIGH    } from '../modules/local/fgbio/filterconsensusreads/main'

include { CALLING_VARDICT                   as CALLINGVARDICTLOW           } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTMED           } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTHIGH          } from '../modules/local/calling_vardict/main'

// include { BBGPOSTANALYSIS                     as BBGPOSTANALYSIS               } from '../modules/local/BBGPOSTANALYSIS/main'

include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTLOW              } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTMED              } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTHIGH             } from '../modules/local/sigprofiler/matrixgenerator/main'

include { SAMTOOLS_FILTER                   as SAMTOOLSFILTERRAW           } from '../modules/local/filter_reads/samtools/main'
include { ASMINUSXS                         as ASMINUSXSRAW                } from '../modules/local/filter_reads/asminusxs/main'
include { SAMTOOLS_FILTER                   as SAMTOOLSFILTERDUPLEX        } from '../modules/local/filter_reads/samtools/main'
include { ASMINUSXS                         as ASMINUSXSDUPLEX             } from '../modules/local/filter_reads/asminusxs/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTP                             as FASTP                       } from '../modules/nf-core/fastp/main'

include { FASTQC                            as PRETRIMFASTQC               } from '../modules/nf-core/fastqc/main'
include { FASTQC                            as FASTQC                      } from '../modules/nf-core/fastqc/main'

include { BAMUTIL_TRIMBAM                   as TRIMBAM                     } from '../modules/nf-core/bamutil/trimbam/main'

include { PICARD_BEDTOINTERVALLIST          as BEDTOINTERVAL               } from '../modules/nf-core/picard/bedtointervallist/main'

//  Metrics
include { QUALIMAP_BAMQC                    as QUALIMAPQC                  } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQC2                 } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCDUPLEX            } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCDUPLEX2           } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCHIGH              } from '../modules/nf-core/qualimap/bamqc/main'

include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHLOW             } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHMED             } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHHIGH            } from '../modules/nf-core/samtools/depth/main'
// include { PICARD_COLLECTMULTIPLEMETRICS     as COLLECTMULTIPLEMETRICS      } from '../modules/nf-core/picard/collectmultiplemetrics/main'

// Versions and reports
include { MULTIQC                                                          } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                      } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Sorting
include { SAMTOOLS_SORT                     as SORTBAM                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCLEAN                } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCONS                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSLOW        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSMED        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSHIGH       } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEX               } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCLEAN          } from '../modules/nf-core/samtools/sort/main'

// include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/nf-core/fgbio/fastqtobam/main'

include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMI             } from '../modules/nf-core/fgbio/groupreadsbyumi/main'
include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { FGBIO_CALLMOLECULARCONSENSUSREADS as CALLMOLECULARCONSENSUSREADS } from '../modules/nf-core/fgbio/callmolecularconsensusreads/main'
include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUSREADS    } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
// include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/nf-core/fgbio/filterconsensusreads/main'
// include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'


// Download annotation cache if needed
include { PREPARE_CACHE                                                   } from '../subworkflows/local/prepare_cache/main'

// Postprocessing of the BAM and the VCF
include { RECOUNT_MUTS                      as RECOUNTMUTSLOW          } from '../subworkflows/local/recount_muts/main'
include { RECOUNT_MUTS                      as RECOUNTMUTSMED          } from '../subworkflows/local/recount_muts/main'
include { RECOUNT_MUTS                      as RECOUNTMUTSHIGH         } from '../subworkflows/local/recount_muts/main'

// Annotation
include { VCF_ANNOTATE_ALL                  as VCFANNOTATELOW          } from '../subworkflows/local/vcf_annotate_all/main'
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEMED          } from '../subworkflows/local/vcf_annotate_all/main'
include { VCF_ANNOTATE_ALL                  as VCFANNOTATEHIGH         } from '../subworkflows/local/vcf_annotate_all/main'










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
        PRETRIMFASTQC(
            INPUT_CHECK.out.reads
        )
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
    // TODO
    // test with real samples whether we could change the "false" here into "true"
    // this would activate sorting the files
    // and would reduce the size of the files stored in the work directory.
    // it works with the test samples
    ALIGNRAWBAM(bam_to_align, ch_ref_index_dir, false)
    ch_versions = ch_versions.mix(ALIGNRAWBAM.out.versions.first())

    SORTBAM(ALIGNRAWBAM.out.bam)
    ch_versions = ch_versions.mix(SORTBAM.out.versions.first())

    // join the bam and the bamindex channels to have
    // the ones from the same samples together
    SORTBAM.out.bam
    .join( SORTBAM.out.csi )
    .set { bam_n_index }

    // QUALIMAPQC(SORTBAM.out.bam, params.targetsfile)

    ASMINUSXSRAW(bam_n_index)
    SAMTOOLSFILTERRAW(ASMINUSXSRAW.out.bam)

    // try to make it template coordinate sorted so that groupbyumi does not need to sort it.
    // we should add the follwoing option in the sam command
    // --template-coordinate
    SORTBAMCLEAN(SAMTOOLSFILTERRAW.out.bam)

    // join the bam and the bamindex channels to have
    // the ones from the same samples together
    SORTBAMCLEAN.out.bam
        .join( SORTBAMCLEAN.out.csi )
        .set { bam_n_index_clean }

    // COLLECTMULTIPLEMETRICS(SORTBAM.out.bam, SORTBAM.out.csi.map{it -> it [1]}, ch_ref_fasta, ch_ref_fasta_fai_index)
    // ch_versions = ch_versions.mix(COLLECTMULTIPLEMETRICS.out.versions.first())

    if (params.targetsfile){
        QUALIMAPQC(SORTBAM.out.bam, params.targetsfile)
        ch_versions = ch_versions.mix(QUALIMAPQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())

        targets_bed = Channel.of([ [ id:"${file(params.targetsfile).getSimpleName()}" ], file(params.targetsfile) ])
        BEDTOINTERVAL(targets_bed, ch_ref_fasta_dict, [])
        ch_versions = ch_versions.mix(BEDTOINTERVAL.out.versions.first())



        // truncate BAM to keep only the reads that are on target
        // TODO
        // see how BAMFILTERREADS requires the BAM file sorted....
        if (params.remove_offtargets){
            BAM_FILTER_READS(bam_n_index_clean,
                            BEDTOINTERVAL.out.interval_list.first().map{it -> it [1]})
            ch_versions = ch_versions.mix(BAM_FILTER_READS.out.versions.first())

            bam_to_group = BAM_FILTER_READS.out.bam
        } else {
            bam_to_group = SORTBAMCLEAN.out.bam
        }


    } else {
        QUALIMAPQC(SORTBAM.out.bam, [])
        ch_versions = ch_versions.mix(QUALIMAPQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())

        bam_to_group = SORTBAMCLEAN.out.bam

    }



    if (params.duplex_seq) {
        //
        // Run fgbio Duplex consensus pipeline
        //

        // MODULE: Run fgbio GroupReadsByUmi
        // requires input template coordinate sorted
        // --template-coordinate
        GROUPREADSBYUMIDUPLEX(bam_to_group, "Paired")
        ch_versions = ch_versions.mix(GROUPREADSBYUMIDUPLEX.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(GROUPREADSBYUMIDUPLEX.out.histogram.map{it[1]}.collect())

        // MODULE: Run fgbio CallDuplexConsensusReads
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)
        ch_versions = ch_versions.mix(CALLDUPLEXCONSENSUSREADS.out.versions.first())

        // MODULE: Run fgbio CollecDuplexSeqMetrics
        COLLECTDUPLEXSEQMETRICS(GROUPREADSBYUMIDUPLEX.out.bam)
        ch_versions = ch_versions.mix(COLLECTDUPLEXSEQMETRICS.out.versions.first())

        // Join groupby stats and duplex seq metrics files from the same samples
        GROUPREADSBYUMIDUPLEX.out.histogram
        .join(COLLECTDUPLEXSEQMETRICS.out.metrics)
        .set {metrics_ch}

        // Plot the family size metrics
        FAMILYMETRICS(metrics_ch)
        ch_versions = ch_versions.mix(FAMILYMETRICS.out.versions.first())
        FAMILYMETRICS.out.log.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)

        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(CALLDUPLEXCONSENSUSREADS.out.bam, ch_ref_index_dir, false)

        SORTBAMDUPLEX(ALIGNDUPLEXCONSENSUSBAM.out.bam)

        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMDUPLEX.out.bam
        .join( SORTBAMDUPLEX.out.csi )
        .set { bam_n_index_duplex }

        // QUALIMAPQCDUPLEX(SORTBAMDUPLEX.out.bam, params.targetsfile)

        ASMINUSXSDUPLEX(bam_n_index_duplex)
        SAMTOOLSFILTERDUPLEX(ASMINUSXSDUPLEX.out.bam)
        
        SORTBAMDUPLEXCLEAN(SAMTOOLSFILTERDUPLEX.out.bam)

        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMDUPLEXCLEAN.out.bam
            .join( SORTBAMDUPLEXCLEAN.out.csi )
            .set { bam_n_index_duplex_clean }

        // requires input coordinate sorted
        QUALIMAPQCDUPLEX(SORTBAMDUPLEXCLEAN.out.bam, params.targetsfile)

        //
        // HIGH CONFIDENCE CALLS
        //

        // MODULE: Run fgbio FilterConsensusReads
        // requires input queryname sorted
        FILTERCONSENSUSREADSHIGH(SAMTOOLSFILTERDUPLEX.out.bam, ch_ref_fasta)
        ch_versions = ch_versions.mix(FILTERCONSENSUSREADSHIGH.out.versions.first())

        // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
        CLIPBAMHIGH(FILTERCONSENSUSREADSHIGH.out.bam, ch_ref_fasta)
        ch_versions = ch_versions.mix(CLIPBAMHIGH.out.versions.first())

        // MODULE: Sort BAM file
        SORTBAMDUPLEXCONSHIGH(CLIPBAMHIGH.out.bam)
        
        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMDUPLEXCONSHIGH.out.bam
        .join( SORTBAMDUPLEXCONSHIGH.out.csi )
        .set { cons_high_bam }


        // Compute depth of the consensus reads aligned to the genome
        COMPUTEDEPTHHIGH(SORTBAMDUPLEXCONSHIGH.out.bam)

        // Quality check
        QUALIMAPQCHIGH(SORTBAMDUPLEXCONSHIGH.out.bam, params.targetsfile)
        ch_versions = ch_versions.mix(QUALIMAPQCHIGH.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCHIGH.out.results.map{it[1]}.collect())


        // Mutation calling for duplex reads
        CALLINGVARDICTHIGH(cons_high_bam,
                            params.targetsfile,
                            ch_ref_fasta, ch_ref_index_dir)
        ch_versions = ch_versions.mix(CALLINGVARDICTHIGH.out.versions.first())

        // Postprocessing the BAM file to get exact coverage per position and allele
        //    also get the Ns per position
        RECOUNTMUTSHIGH(cons_high_bam,
                        CALLINGVARDICTHIGH.out.vcf,
                        params.targetsfile,
                        ch_ref_fasta
                        )
        ch_versions = ch_versions.mix(RECOUNTMUTSHIGH.out.versions.first())


        VCFANNOTATEHIGH(CALLINGVARDICTHIGH.out.vcf,
                            ch_ref_fasta,
                            params.vep_genome,
                            params.vep_species,
                            params.vep_cache_version,
                            vep_cache,
                            vep_extra_files)
        ch_versions = ch_versions.mix(VCFANNOTATEHIGH.out.versions.first())

        CALLINGVARDICTHIGH.out.vcf.map{it -> it[1]}.set { mutation_files_high }

        SIGPROFPLOTHIGH(mutation_files_high.collect())
        // ch_versions = ch_versions.mix(SIGPROFPLOTDUPLEX.out.versions.first())


        // TODO
        // add a module that plots the signature, with the correct normalization
        // use bgsignature, and the input targets file

        //
        // MEDIUM CONFIDENCE CALLS
        //
        if (params.duplex_med_conf){

            FILTERCONSENSUSREADSMED(SAMTOOLSFILTERDUPLEX.out.bam, ch_ref_fasta)

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMMED(FILTERCONSENSUSREADSMED.out.bam, ch_ref_fasta)
            ch_versions = ch_versions.mix(CLIPBAMMED.out.versions.first())

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONSMED(CLIPBAMMED.out.bam)

            // join the bam and the bamindex channels to have
            // the ones from the same samples together
            SORTBAMDUPLEXCONSMED.out.bam
            .join( SORTBAMDUPLEXCONSMED.out.csi )
            .set { cons_med_bam }

            // Compute depth of the consensus reads aligned to the genome
            COMPUTEDEPTHMED(SORTBAMDUPLEXCONSMED.out.bam)

            // Mutation calling for all reads
            CALLINGVARDICTMED(cons_med_bam,
                                params.targetsfile,
                                ch_ref_fasta, ch_ref_index_dir)

            // Postprocessing the BAM file to get exact coverage per position and allele
            //    also get the Ns per position
            RECOUNTMUTSMED(cons_med_bam,
                            CALLINGVARDICTMED.out.vcf,
                            params.targetsfile,
                            ch_ref_fasta,
                            params.filter_mutations // another option would be to pass this as a tuple of as many elements
                                                //      as possible filters so that we can tune which filters are applied
                                                )
            ch_versions = ch_versions.mix(RECOUNTMUTSMED.out.versions.first())

            VCFANNOTATEMED(CALLINGVARDICTMED.out.vcf,
                            ch_ref_fasta,
                            params.vep_genome,
                            params.vep_species,
                            params.vep_cache_version,
                            vep_cache,
                            vep_extra_files)

            CALLINGVARDICTMED.out.vcf.map{it -> it[1]}.set { mutation_files_med }
            SIGPROFPLOTMED(mutation_files_med.collect())
        }



        //
        // Low Confidence mutations
        //
        if (params.duplex_low_conf){

            // filter the reads with the low conf parameters
            FILTERCONSENSUSREADSLOW(SAMTOOLSFILTERDUPLEX.out.bam, ch_ref_fasta)

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMLOW(FILTERCONSENSUSREADSLOW.out.bam, ch_ref_fasta)
            ch_versions = ch_versions.mix(CLIPBAMLOW.out.versions.first())

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONSLOW(CLIPBAMLOW.out.bam)

            // join the bam and the bamindex channels to have
            // the ones from the same samples together
            SORTBAMDUPLEXCONSLOW.out.bam
            .join( SORTBAMDUPLEXCONSLOW.out.csi )
            .set { cons_low_bam }

            // Compute depth of the consensus reads aligned to the genome
            COMPUTEDEPTHLOW(SORTBAMDUPLEXCONSLOW.out.bam)

            // Mutation calling for all reads
            CALLINGVARDICTLOW(cons_low_bam,
                                params.targetsfile,
                                ch_ref_fasta, ch_ref_index_dir)

            // Postprocessing the BAM file to get exact coverage per position and allele
            //    also get the Ns per position
            RECOUNTMUTSLOW(cons_low_bam,
                            CALLINGVARDICTLOW.out.vcf,
                            params.targetsfile,
                            ch_ref_fasta,
                            params.filter_mutations // another option would be to pass this as a tuple of as many elements
                                                //      as possible filters so that we can tune which filters are applied
                                                )
            ch_versions = ch_versions.mix(RECOUNTMUTSMED.out.versions.first())

            VCFANNOTATELOW(CALLINGVARDICTLOW.out.vcf,
                            ch_ref_fasta,
                            params.vep_genome,
                            params.vep_species,
                            params.vep_cache_version,
                            vep_cache,
                            vep_extra_files)

            CALLINGVARDICTLOW.out.vcf.map{it -> it[1]}.set { mutation_files_low }
            SIGPROFPLOTLOW(mutation_files_low.collect())
        }



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
        SORTBAMCONS.out.bam
        .join(SORTBAMCONS.out.csi)
        .set {umi_bam}

        // Mutation calling for non-duplex reads
        CALLINGVARDICT(umi_bam,
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
