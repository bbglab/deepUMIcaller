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

include { SAMTOOLS_FILTER                   as SAMTOOLSFILTERDUPLEX             } from '../modules/local/filter_reads/samtools/main'
include { ASMINUSXS                         as ASMINUSXSDUPLEX                  } from '../modules/local/filter_reads/asminusxs/main'

include { FGBIO_CLIPBAM                     as CLIPBAMLOW                       } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMMED                       } from '../modules/local/clipbam/main'
include { FGBIO_CLIPBAM                     as CLIPBAMHIGH                      } from '../modules/local/clipbam/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSAM           } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSLOW          } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSMED          } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSHIGH         } from '../modules/local/fgbio/filterconsensusreads/main'

include { CREATEBED_FROM_TSV                as CREATEBEDLOW                     } from '../modules/local/createbed/main'
include { CREATEBED_FROM_TSV                as CREATEBEDMED                     } from '../modules/local/createbed/main'
include { CREATEBED_FROM_TSV                as CREATEBEDHIGH                    } from '../modules/local/createbed/main'

include { CALLING_VARDICT                   as CALLINGVARDICTLOW                } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTMED                } from '../modules/local/calling_vardict/main'
include { CALLING_VARDICT                   as CALLINGVARDICTHIGH               } from '../modules/local/calling_vardict/main'

include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTLOW                   } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTLOWPUR               } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTLOWPYR               } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTMED                   } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTMEDPUR               } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTMEDPYR               } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTHIGH                  } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTHIGHPUR               } from '../modules/local/sigprofiler/matrixgenerator/main'
include { SIGPROFILER_MATRIXGENERATOR       as SIGPROFPLOTHIGHPYR               } from '../modules/local/sigprofiler/matrixgenerator/main'


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
include { QUALIMAP_BAMQC                    as QUALIMAPQCDUPLEX            } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCHIGH              } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCMED               } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC                    as QUALIMAPQCLOW               } from '../modules/nf-core/qualimap/bamqc/main'


include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHLOW             } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHMED             } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_DEPTH                    as COMPUTEDEPTHHIGH            } from '../modules/nf-core/samtools/depth/main'

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
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSLOW        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSMED        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONSHIGH       } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEX               } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCLEAN          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXFILTERED       } from '../modules/nf-core/samtools/sort/main'

// include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                  } from '../modules/nf-core/fgbio/fastqtobam/main'

include { FGBIO_GROUPREADSBYUMI             as GROUPREADSBYUMIDUPLEX       } from '../modules/nf-core/fgbio/groupreadsbyumi/main'

include { FGBIO_CALLDUPLEXCONSENSUSREADS    as CALLDUPLEXCONSENSUSREADS    } from '../modules/nf-core/fgbio/callduplexconsensusreads/main'
// include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADS        } from '../modules/nf-core/fgbio/filterconsensusreads/main'
// include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS     } from '../modules/nf-core/fgbio/collectduplexseqmetrics/main'


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

workflow DEEPUMICALLER {

    
    if (params.ref_fasta) {
        ch_ref_fasta = Channel.fromPath(params.ref_fasta).collect()
    
        // define additional fasta file names
        ch_ref_fasta_file = file(params.ref_fasta, checkIfExists: true)
        ch_ref_fasta_dict = file("${ch_ref_fasta_file.parent/ch_ref_fasta_file.baseName}.dict", checkIfExists: true)
    } else {
        log.error "No reference FASTA was specified (--ref_fasta)."
        exit 1
    }
    
    ch_ref_index_dir = ch_ref_fasta.map { it -> it.parent }

    ch_multiqc_files = Channel.empty()


    vep_cache = params.vep_cache
    vep_extra_files = []
    

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
                QUALIMAPQC(SORTBAM.out.bam, params.targetsfile)
                
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())
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
                                BEDTOINTERVAL.out.interval_list.first().map{it -> it [1]})
                

                bam_to_group = BAM_FILTER_READS.out.bam
            } else {
                bam_to_group = SORTBAMCLEAN.out.bam
            }


        } else {
            if (params.perform_qcs){
                QUALIMAPQC(SORTBAM.out.bam, [])
                
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQC.out.results.map{it[1]}.collect())
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
        
        FAMILYMETRICS.out.sample_data.map{it -> it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)
        FAMILYMETRICS.out.curve_data.map{it -> it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CollecDuplexSeqMetrics only on target
        COLLECTDUPLEXSEQMETRICSONTARGET(GROUPREADSBYUMIDUPLEX.out.bam, BEDTOINTERVAL.out.interval_list.first().map{it -> it[1]} )
        


        // Plot the family size metrics
        FAMILYMETRICSONTARGET(COLLECTDUPLEXSEQMETRICSONTARGET.out.metrics)
        
        FAMILYMETRICSONTARGET.out.sample_data.map{it -> it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)
        FAMILYMETRICSONTARGET.out.curve_data.map{it -> it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CallDuplexConsensusReads
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)
        


        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(CALLDUPLEXCONSENSUSREADS.out.bam, ch_ref_index_dir, false)


        SORTBAMDUPLEX(ALIGNDUPLEXCONSENSUSBAM.out.bam)

        // join the bam and the bamindex channels to have
        // the ones from the same samples together
        SORTBAMDUPLEX.out.bam
        .join( SORTBAMDUPLEX.out.csi )
        .set { bam_n_index_duplex }

        ASMINUSXSDUPLEX(bam_n_index_duplex)
        SAMTOOLSFILTERDUPLEX(ASMINUSXSDUPLEX.out.bam)
        SORTBAMDUPLEXFILTERED(SAMTOOLSFILTERDUPLEX.out.bam)

        duplex_filtered_bam = SORTBAMDUPLEXFILTERED.out.bam

        ASMINUSXSDUPLEX.out.discarded_bam.map{it -> [it[0], params.targetsfile, it[1]]}.set { discarded_bam_targeted }
        DISCARDEDCOVERAGETARGETED(discarded_bam_targeted, [])
        

        ASMINUSXSDUPLEX.out.discarded_bam.map{it -> [it[0], params.global_exons_file, it[1]]}.set { discarded_bam }
        DISCARDEDCOVERAGEGLOBAL(discarded_bam, [])
        

    }

    if (params.step == 'filterconsensus') {
        duplex_filtered_bam = INPUT_CHECK.out.reads
    }


    FILTERCONSENSUSREADSAM(duplex_filtered_bam, ch_ref_fasta)
    SORTBAMDUPLEXCLEAN(FILTERCONSENSUSREADSAM.out.bam)
    
    // join the bam and the bamindex channels to have
    // the ones from the same samples together
    SORTBAMDUPLEXCLEAN.out.bam
    .join( SORTBAMDUPLEXCLEAN.out.csi )
    .set { bam_n_index_duplex_clean }

    if (params.perform_qcs){
        // requires input coordinate sorted
        QUALIMAPQCDUPLEX(SORTBAMDUPLEXCLEAN.out.bam, params.targetsfile)
        ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCDUPLEX.out.results.map{it[1]}.collect())
    }

    //
    // HIGH CONFIDENCE CALLS
    //
    if (params.duplex_high_conf) {

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus']) {

            // MODULE: Run fgbio FilterConsensusReads
            // requires input queryname sorted
            FILTERCONSENSUSREADSHIGH(duplex_filtered_bam, ch_ref_fasta)
            

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMHIGH(FILTERCONSENSUSREADSHIGH.out.bam, ch_ref_fasta)
            

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONSHIGH(CLIPBAMHIGH.out.bam)
            

            // join the bam and the bamindex channels to have
            // the ones from the same samples together
            SORTBAMDUPLEXCONSHIGH.out.bam
            .join( SORTBAMDUPLEXCONSHIGH.out.csi )
            .set { cons_high_bam }

            if (params.perform_qcs){
                // Quality check
                QUALIMAPQCHIGH(SORTBAMDUPLEXCONSHIGH.out.bam, params.targetsfile)
                
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCHIGH.out.results.map{it[1]}.collect())
            }
        }

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus', 'calling']) {

            if (params.step == 'calling') {
                cons_high_bam = INPUT_CHECK.out.reads
            }

            cons_high_bam.map{it -> [it[0], it[1]] }
            .set{ cons_high_bam_only }

            // Compute depth of the consensus reads aligned to the genome
            COMPUTEDEPTHHIGH(cons_high_bam_only)
            

            CREATEBEDHIGH(COMPUTEDEPTHHIGH.out.tsv)
            

            cons_high_bam
            .join( CREATEBEDHIGH.out.bed )
            .set { cons_high_bam_bed }

            // Mutation calling for duplex reads
            CALLINGVARDICTHIGH(cons_high_bam_bed,
                                ch_ref_fasta, ch_ref_index_dir)
            

            // Postprocessing the BAM file to get exact coverage per position and allele
            //    also get the Ns per position
            RECOUNTMUTSHIGH(cons_high_bam,
                            bam_n_index_duplex_clean,
                            CALLINGVARDICTHIGH.out.vcf,
                            CREATEBEDHIGH.out.bed,
                            ch_ref_fasta
                        )
            


            if (params.annotate_mutations){
                VCFANNOTATEHIGH(CALLINGVARDICTHIGH.out.vcf,
                                ch_ref_fasta,
                                vep_cache,
                                vep_extra_files)
                
            }

            RECOUNTMUTSHIGH.out.somatic_vcf.map{it -> it[1]}.set { mutation_files_high }
            SIGPROFPLOTHIGH(mutation_files_high.collect())
            

            RECOUNTMUTSHIGH.out.purvcf.map{it -> it[1]}.set { mutation_files_pur_high }
            SIGPROFPLOTHIGHPUR(mutation_files_pur_high.collect())
            

            RECOUNTMUTSHIGH.out.pyrvcf.map{it -> it[1]}.set { mutation_files_pyr_high }
            SIGPROFPLOTHIGHPYR(mutation_files_pyr_high.collect())
            

        }

    }


    //
    // MEDIUM CONFIDENCE CALLS
    //
    if (params.duplex_med_conf){

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus']) {
            FILTERCONSENSUSREADSMED(duplex_filtered_bam, ch_ref_fasta)

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMMED(FILTERCONSENSUSREADSMED.out.bam, ch_ref_fasta)
            

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONSMED(CLIPBAMMED.out.bam)

            // join the bam and the bamindex channels to have
            // the ones from the same samples together
            SORTBAMDUPLEXCONSMED.out.bam
            .join( SORTBAMDUPLEXCONSMED.out.csi )
            .set { cons_med_bam }

            // Quality check
            if (params.perform_qcs){
                QUALIMAPQCMED(SORTBAMDUPLEXCONSMED.out.bam, params.targetsfile)
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCMED.out.results.map{it[1]}.collect())
                
                SORTBAMDUPLEXCONSMED.out.bam.map{it -> [it[0], params.global_exons_file, it[1]]}.set { duplex_filt_bam_n_bed }
                COVERAGEGLOBAL(duplex_filt_bam_n_bed, [])

            }
        }

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus', 'calling']) {

            // ASSIGN cons_med_bam = to our input bam
            if (params.step == 'calling') {
                cons_med_bam = INPUT_CHECK.out.reads
            }

            cons_med_bam.map{it -> [it[0], it[1]] }
            .set{ cons_med_bam_only }

            // Compute depth of the consensus reads aligned to the genome
            COMPUTEDEPTHMED(cons_med_bam_only)

            CREATEBEDMED(COMPUTEDEPTHMED.out.tsv)

            cons_med_bam
            .join( CREATEBEDMED.out.bed )
            .set { cons_med_bam_bed }

            // Mutation calling for all reads
            CALLINGVARDICTMED(cons_med_bam_bed,
                                ch_ref_fasta, ch_ref_index_dir)

            // Postprocessing the BAM file to get exact coverage per position and allele
            //    also get the Ns per position
            RECOUNTMUTSMED(cons_med_bam,
                            bam_n_index_duplex_clean,
                            CALLINGVARDICTMED.out.vcf,
                            CREATEBEDMED.out.bed,
                            ch_ref_fasta)
            

            if (params.annotate_mutations){
                VCFANNOTATEMED(CALLINGVARDICTMED.out.vcf,
                                ch_ref_fasta,
                                vep_cache,
                                vep_extra_files)
            }

            RECOUNTMUTSMED.out.somatic_vcf.map{it -> it[1]}.set { mutation_files_med }
            SIGPROFPLOTMED(mutation_files_med.collect())

            RECOUNTMUTSMED.out.purvcf.map{it -> it[1]}.set { mutation_files_pur_med }
            SIGPROFPLOTMEDPUR(mutation_files_pur_med.collect())
            

            RECOUNTMUTSMED.out.pyrvcf.map{it -> it[1]}.set { mutation_files_pyr_med }
            SIGPROFPLOTMEDPYR(mutation_files_pyr_med.collect())
            

        }
    }

    //
    // Low Confidence mutations
    //
    if (params.duplex_low_conf){

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus']) {


            // filter the reads with the low conf parameters
            FILTERCONSENSUSREADSLOW(duplex_filtered_bam, ch_ref_fasta)

            // MODULE: Hard clipping read pairs that overlap, and that go beyond the pair starting point
            CLIPBAMLOW(FILTERCONSENSUSREADSLOW.out.bam, ch_ref_fasta)
            

            // MODULE: Sort BAM file
            SORTBAMDUPLEXCONSLOW(CLIPBAMLOW.out.bam)

            // join the bam and the bamindex channels to have
            // the ones from the same samples together
            SORTBAMDUPLEXCONSLOW.out.bam
            .join( SORTBAMDUPLEXCONSLOW.out.csi )
            .set { cons_low_bam }
            
            // Quality check
            if (params.perform_qcs){
                QUALIMAPQCLOW(SORTBAMDUPLEXCONSLOW.out.bam, params.targetsfile)
                
                ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCLOW.out.results.map{it[1]}.collect())
            }
        }

        if (params.step in ['mapping', 'groupreadsbyumi', 'filterconsensus', 'calling']) {

            if (params.step == 'calling') {
                cons_low_bam = INPUT_CHECK.out.reads
            }

            cons_low_bam.map{it -> [it[0], it[1]] }
            .set{ cons_low_bam_only }

            // Compute depth of the consensus reads aligned to the genome
            COMPUTEDEPTHLOW(cons_low_bam_only)

            CREATEBEDLOW(COMPUTEDEPTHLOW.out.tsv)

            cons_low_bam
            .join( CREATEBEDLOW.out.bed )
            .set { cons_low_bam_bed }


            // Mutation calling for all reads
            CALLINGVARDICTLOW(cons_low_bam_bed,
                                ch_ref_fasta, ch_ref_index_dir)

            // Postprocessing the BAM file to get exact coverage per position and allele
            //    also get the Ns per position
            RECOUNTMUTSLOW(cons_low_bam,
                            bam_n_index_duplex_clean,
                            CALLINGVARDICTLOW.out.vcf,
                            CREATEBEDLOW.out.bed,
                            ch_ref_fasta)
            

            if (params.annotate_mutations){
                VCFANNOTATELOW(CALLINGVARDICTLOW.out.vcf,
                                ch_ref_fasta,
                                vep_cache,
                                vep_extra_files)
            }

            RECOUNTMUTSLOW.out.somatic_vcf.map{it -> it[1]}.set { mutation_files_low }
            SIGPROFPLOTLOW(mutation_files_low.collect())

            RECOUNTMUTSLOW.out.purvcf.map{it -> it[1]}.set { mutation_files_pur_low }
            SIGPROFPLOTLOWPUR(mutation_files_pur_low.collect())

            RECOUNTMUTSLOW.out.pyrvcf.map{it -> it[1]}.set { mutation_files_pyr_low }
            SIGPROFPLOTLOWPYR(mutation_files_pyr_low.collect())
            

        }

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
