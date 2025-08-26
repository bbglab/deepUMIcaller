/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK                                                           } from '../subworkflows/local/input_check'

include { BAM_FILTER_READS                                                      } from '../subworkflows/local/bam_filter_reads/main'


include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                       } from '../modules/local/fgbio/fastqtobam/main'

include { ALIGN_BAM                         as ALIGNRAWBAM                      } from '../modules/local/align_bam/main'

include { ALIGN_BAM                         as ALIGNCONSENSUSBAM                } from '../modules/local/align_bam/main'
include { ALIGN_BAM                         as ALIGNDUPLEXCONSENSUSBAM          } from '../modules/local/align_bam/main'

include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS          } from '../modules/local/fgbio/collectduplexseqmetrics/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICSONTARGET  } from '../modules/local/fgbio/collectduplexseqmetrics/main'

include { FAMILYSIZEMETRICS                 as FAMILYMETRICS                    } from '../modules/local/familymetrics/main'
include { FAMILYSIZEMETRICS                 as FAMILYMETRICSONTARGET            } from '../modules/local/familymetrics/main'

include { SAMTOOLS_FILTER                   as SAMTOOLSFILTERALLMOLECULES       } from '../modules/local/filter_reads/samtools/main'
include { ASMINUSXS                         as ASMINUSXSDUPLEX                  } from '../modules/local/filter_reads/asminusxs/main'

include { FGBIO_CLIPBAM                     as CLIPBAM                          } from '../modules/local/clipbam/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSAM           } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSDUPLEX       } from '../modules/local/fgbio/filterconsensusreads/main'

include { CREATEBED_FROM_TSV                as CREATEBED                        } from '../modules/local/createbed/main'

include { CALLING_VARDICT                   as CALLINGVARDICTDUPLEX             } from '../modules/local/calling_vardict/main'

include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOT                      } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTPUR                   } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTPYR                   } from '../modules/local/sigprofiler/matrixgenerator/main'


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
include { QUALIMAP_BAMQC                    as QUALIMAPQCRAW               } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCALLMOLECULES      } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCDUPLEX            } from '../modules/nf-core/qualimap/bamqc/main'


include { SAMTOOLS_DEPTH                    as COMPUTEDEPTH                } from '../modules/nf-core/samtools/depth/main'

include { BEDTOOLS_COVERAGE                 as DISCARDEDCOVERAGETARGETED   } from '../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_COVERAGE                 as DISCARDEDCOVERAGEGLOBAL     } from '../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_COVERAGE                 as COVERAGEGLOBAL              } from '../modules/nf-core/bedtools/coverage/main'


// Versions and reports
include { MULTIQC                                                          } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                      } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Sorting
include { SAMTOOLS_SORT                     as SORTBAM                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCLEAN                } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCONS                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONS           } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMALLMOLECULES         } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMAMCLEAN              } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMAMFILTERED           } from '../modules/nf-core/samtools/sort/main'

// include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/nf-core/fgbio/fastqtobam/main'

include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUSREADS    } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
// include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/nf-core/fgbio/filterconsensusreads/main'
// include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'


// Postprocessing of the BAM and the VCF
include { RECOUNT_MUTS                      as RECOUNTMUTS                 } from '../subworkflows/local/recount_muts/main'

// Annotation
include { VCF_ANNOTATE_ALL                  as VCFANNOTATE                 } from '../subworkflows/local/vcf_annotate_all/main'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DEEPUMICALLER {

    
    if (params.ref_fasta) {
        ch_ref_fasta = file(params.ref_fasta, checkIfExists: true)
        ch_ref_fasta_dict = file("${ch_ref_fasta.parent/ch_ref_fasta.baseName}.dict", checkIfExists: true)
        ch_ref_index_dir = ch_ref_fasta.parent
    } else {
        log.error "No reference FASTA was specified (--ref_fasta)."
        exit 1
    }

    ch_multiqc_files = Channel.empty()
   

    if (params.targetsfile) {
        targets_bed = Channel.of([ [ id:"${file(params.targetsfile).getSimpleName()}" ], file(params.targetsfile) ])
        BEDTOINTERVAL(targets_bed, ch_ref_fasta_dict, [])
    }


    
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    INPUT_CHECK (
        file(params.input), 
        params.step
    )
    

    if (params.step == 'mapping') {

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
            
            reads_to_qc = FASTP.out.reads
        } else {
            reads_to_qc = INPUT_CHECK.out.reads
        }


        // MODULE: Run FastQC
        FASTQC (
            reads_to_qc
        )
        
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))


        // MODULE: Run fgbio FastqToBam
        // to get the UMIs out of the reads and into the tag
        FASTQTOBAM(reads_to_qc)
        


        // Decide whether we clip the beginning and/or end of the reads or nothing
        if ( (params.left_clip > 0) || (params.right_clip > 0) ) {
            // params.left_clip = 4     ->      Remove 4bp from the 5' end of the reads
            TRIMBAM(FASTQTOBAM.out.bam, params.left_clip, params.right_clip)
            
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
        

        if (params.perform_qcs){
            SORTBAM(ALIGNRAWBAM.out.bam)
            
        }

        // template coordinate sorting for the GroupByUMI
        SORTBAMCLEAN(ALIGNRAWBAM.out.bam)
        


        if (params.targetsfile){
            if (params.perform_qcs){
                QUALIMAPQCRAW(SORTBAM.out.bam, params.targetsfile)
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCRAW.out.results.map{it[1]}.collect())
            }
            // truncate BAM to keep only the reads that are on target
            // TODO
            // see how BAMFILTERREADS requires the BAM file sorted....
            if (params.remove_offtargets){
                // join the bam and the bamindex channels to have
                // the ones from the same samples together
                SORTBAMCLEAN.out.bam
                .join( SORTBAMCLEAN.out.csi )
                .set { bam_n_index_clean }

                BAM_FILTER_READS(bam_n_index_clean,
                                BEDTOINTERVAL.out.interval_list.first().map{it[1]})
                

                bam_to_group = BAM_FILTER_READS.out.bam
            } else {
                bam_to_group = SORTBAMCLEAN.out.bam
            }


        } else {
            if (params.perform_qcs){
                QUALIMAPQCRAW(SORTBAM.out.bam, [])
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCRAW.out.results.map{it[1]}.collect())
            }
            bam_to_group = SORTBAMCLEAN.out.bam
        }
    }


    //
    // Run fgbio Duplex consensus pipeline
    //

    if (params.step in ['mapping', 'groupreadsbyumi']) {

        // ASSIGN bam_to_group = to our input bam
        if (params.step == 'groupreadsbyumi') {
            bam_to_group = INPUT_CHECK.out.reads
        }

        // MODULE: Run fgbio GroupReadsByUmi
        // requires input template coordinate sorted
        GROUPREADSBYUMIDUPLEX(bam_to_group, "Paired")
        
        ch_multiqc_files = ch_multiqc_files.mix(GROUPREADSBYUMIDUPLEX.out.histogram.map{it[1]}.collect())


        // MODULE: Run fgbio CollecDuplexSeqMetrics
        COLLECTDUPLEXSEQMETRICS(GROUPREADSBYUMIDUPLEX.out.bam, [])
        


        // Plot the family size metrics
        FAMILYMETRICS(COLLECTDUPLEXSEQMETRICS.out.metrics)

        FAMILYMETRICS.out.sample_data.map{it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)
        FAMILYMETRICS.out.curve_data.map{it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CollecDuplexSeqMetrics only on target
        COLLECTDUPLEXSEQMETRICSONTARGET(GROUPREADSBYUMIDUPLEX.out.bam, BEDTOINTERVAL.out.interval_list.first().map{it[1]} )
        


        // Plot the family size metrics
        FAMILYMETRICSONTARGET(COLLECTDUPLEXSEQMETRICSONTARGET.out.metrics)

        FAMILYMETRICSONTARGET.out.sample_data.map{it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)
        FAMILYMETRICSONTARGET.out.curve_data.map{it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CallDuplexConsensusReads
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)
        

        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(CALLDUPLEXCONSENSUSREADS.out.bam, ch_ref_index_dir, false)


        SORTBAMALLMOLECULES(ALIGNDUPLEXCONSENSUSBAM.out.bam)

        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMALLMOLECULES.out.bam
        .join( SORTBAMALLMOLECULES.out.csi )
        .set { bam_n_index_all_molecules }

        ASMINUSXSDUPLEX(bam_n_index_all_molecules)
        SAMTOOLSFILTERALLMOLECULES(ASMINUSXSDUPLEX.out.bam)
        SORTBAMAMFILTERED(SAMTOOLSFILTERALLMOLECULES.out.bam)

        duplex_filtered_bam = SORTBAMAMFILTERED.out.bam

        ASMINUSXSDUPLEX.out.discarded_bam.map{[it[0], params.targetsfile, it[1]]}.set { discarded_bam_targeted }
        DISCARDEDCOVERAGETARGETED(discarded_bam_targeted, [])

        ASMINUSXSDUPLEX.out.discarded_bam.map{[it[0], params.global_exons_file, it[1]]}.set { discarded_bam }
        DISCARDEDCOVERAGEGLOBAL(discarded_bam, [])

    }

    if (params.step == 'filterconsensus') {
        duplex_filtered_bam = INPUT_CHECK.out.reads
    }


    FILTERCONSENSUSREADSAM(duplex_filtered_bam, ch_ref_fasta)
    SORTBAMAMCLEAN(FILTERCONSENSUSREADSAM.out.bam)
    
    // join the bam and the bamindex channels to have
    // the ones from the same samples together
    SORTBAMAMCLEAN.out.bam
    .join( SORTBAMAMCLEAN.out.csi )
    .set { bam_n_index_duplex_clean }

    if (params.perform_qcs){
        // requires input coordinate sorted
        QUALIMAPQCALLMOLECULES(SORTBAMAMCLEAN.out.bam, params.targetsfile)
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCALLMOLECULES.out.results.map{it[1]}.collect())
    }

    //
    // DUPLEX CALLS
    //
    if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus']) {
        FILTERCONSENSUSREADSDUPLEX(duplex_filtered_bam, ch_ref_fasta)

        // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
        CLIPBAM(FILTERCONSENSUSREADSDUPLEX.out.bam, ch_ref_fasta)
        

        // MODULE: Sort BAM file
        SORTBAMDUPLEXCONS(CLIPBAM.out.bam)

        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMDUPLEXCONS.out.bam
        .join( SORTBAMDUPLEXCONS.out.csi )
        .set { cons_duplex_bam }

        // Quality check
        if (params.perform_qcs){
            QUALIMAPQCDUPLEX(SORTBAMDUPLEXCONS.out.bam, params.targetsfile)
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCDUPLEX.out.results.map{it[1]}.collect())

            SORTBAMDUPLEXCONS.out.bam.map{it -> [it[0], params.global_exons_file, it[1]]}.set { duplex_filt_bam_n_bed }
            COVERAGEGLOBAL(duplex_filt_bam_n_bed, [])
        }
    }

    if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus', 'calling']) {
    
        // ASSIGN cons_med_bam = to our input bam
        if (params.step == 'calling') {
            cons_duplex_bam = INPUT_CHECK.out.reads
        }

        cons_duplex_bam.map{[it[0], it[1]]}
        .set{ cons_duplex_bam_only }

        // Compute depth of the consensus reads aligned to the genome
        COMPUTEDEPTH(cons_duplex_bam_only)

        CREATEBED(COMPUTEDEPTH.out.tsv)

        cons_duplex_bam
        .join( CREATEBED.out.bed )
        .set { cons_duplex_bam_bed }

        // Mutation calling for all reads
        CALLINGVARDICTDUPLEX(cons_duplex_bam_bed,
                                ch_ref_fasta, ch_ref_index_dir)

        // Postprocessing the BAM file to get exact coverage per position and allele
        //    also get the Ns per position
        RECOUNTMUTS(cons_duplex_bam,
                        bam_n_index_duplex_clean,
                        CALLINGVARDICTDUPLEX.out.vcf,
                        CREATEBED.out.bed,
                        ch_ref_fasta)

        if (params.annotate_mutations){
            VCFANNOTATE(CALLINGVARDICTDUPLEX.out.vcf,
                            ch_ref_fasta)
        }

        RECOUNTMUTS.out.somatic_vcf.map{it[1]}.set { mutation_files_duplex }
        SIGPROFPLOT(mutation_files_duplex.collect())

        RECOUNTMUTS.out.purvcf.map{it[1]}.set { mutation_files_pur_duplex }
        SIGPROFPLOTPUR(mutation_files_pur_duplex.collect())

        RECOUNTMUTS.out.pyrvcf.map{it[1]}.set { mutation_files_pyr_duplex }
        SIGPROFPLOTPYR(mutation_files_pyr_duplex.collect())

    }


    CUSTOM_DUMPSOFTWAREVERSIONS (
        Channel.topic('versions').unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    def summary_params  = paramsSummaryMap(workflow)
    workflow_summary    = paramsSummaryMultiqc(summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)


    ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
    
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect()
    )
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
