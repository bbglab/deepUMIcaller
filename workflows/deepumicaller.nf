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

include { SPLITFASTQ                                                            } from '../modules/local/splitfastq/main'

include { FGBIO_FASTQTOBAM                  as FASTQTOBAM                       } from '../modules/local/fgbio/fastqtobam/main'

include { ALIGN_BAM                         as ALIGNRAWBAM                      } from '../modules/local/align_bam/main'

include { MERGEBAM                                                              } from '../modules/local/mergebam/main'
include { MERGEBAM                          as MERGEBAMDUPLEX                   } from '../modules/local/mergebam/main'
include { SPLITBAMCHROM                                                         } from '../modules/local/splitbamchrom/main'

include { ALIGN_BAM                         as ALIGNCONSENSUSBAM                } from '../modules/local/align_bam/main'
include { ALIGN_BAM                         as ALIGNDUPLEXCONSENSUSBAM          } from '../modules/local/align_bam/main'

include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICS          } from '../modules/local/fgbio/collectduplexseqmetrics/main'
include { FGBIO_COLLECTDUPLEXSEQMETRICS     as COLLECTDUPLEXSEQMETRICSONTARGET  } from '../modules/local/fgbio/collectduplexseqmetrics/main'

include { FAMILYSIZEMETRICS                 as FAMILYMETRICS                    } from '../modules/local/familymetrics/main'
include { FAMILYSIZEMETRICS                 as FAMILYMETRICSONTARGET            } from '../modules/local/familymetrics/main'

include { UNMAP_BAM                         as UNMAPBAM                         } from '../modules/local/unmap/main'
include { SAMTOOLS_FILTER                   as SAMTOOLSFILTERALLMOLECULES       } from '../modules/local/filter_reads/samtools/main'
include { ASMINUSXS                         as ASMINUSXSDUPLEX                  } from '../modules/local/filter_reads/asminusxs/main'

include { FGBIO_CLIPBAM                     as CLIPBAM                          } from '../modules/local/clipbam/main'

include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSAM           } from '../modules/local/fgbio/filterconsensusreads/main'
include { FGBIO_FILTERCONSENSUSREADS        as FILTERCONSENSUSREADSDUPLEX       } from '../modules/local/fgbio/filterconsensusreads/main'

include { CREATEBED_FROM_TSV                as CREATEBED                        } from '../modules/local/createbed/main'

include { BAM_CALL_VARDICT_PARALLEL                                             } from '../subworkflows/local/bam_call_vardict_parallel/main'

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

include { PICARD_MERGESAMFILES              as MERGEBAMS                    } from '../modules/nf-core/picard/mergesamfiles/main'

// Versions and reports
include { MULTIQC                                                          } from '../modules/nf-core/multiqc/main'
include { MULTIQC                           as MULTIQCDUPLEX               } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                                      } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// Sorting
include { SAMTOOLS_SORT                     as SORTBAM                     } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCLEAN                } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMCONS                 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMDUPLEXCONS           } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMALLMOLECULES         } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMAMCLEAN              } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMAMFILTERED           } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT                     as SORTBAMMERGED               } from '../modules/nf-core/samtools/sort/main'

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
    ch_multiqc_config = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
   

    // Create value channels for targets and global exons files (if provided)
    ch_targetsfile = params.targetsfile ? file(params.targetsfile, checkIfExists: true) : Channel.empty()
    ch_global_exons_file = params.global_exons_file ? file(params.global_exons_file, checkIfExists: true) : file(params.targetsfile, checkIfExists: true)


    targets_bed = Channel.of([ [ id:"${file(params.targetsfile).getSimpleName()}" ], ch_targetsfile ])
    BEDTOINTERVAL(targets_bed, ch_ref_fasta_dict, [])


    
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

        // Optional: Include fastqs split
        if (params.run_splitfastq) {
            SPLITFASTQ ( reads_to_qc )
            SPLITFASTQ.out.split_fastqs
            .transpose()
            .map { meta, fastqs -> 
                def new_meta = meta.clone()
                // Extract part ID from FASTQ filename using regex, e.g. "sample1_L001_R1_001.fastq.gz" -> "sample1"
                def match = fastqs[0].name =~ /^([^. _-]+)/
                def part_id = match ? match[0][1] : "unknown"
                new_meta.id = "${meta.id}_${part_id}"
                [new_meta, fastqs]
            }
            .set { split_fastqs_ch }
        } else {
            split_fastqs_ch = reads_to_qc
        }

        FASTQTOBAM(split_fastqs_ch)


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
        //ALIGNRAWBAM.out.bam.view { meta, bam -> "Debug ALIGNRAWBAM: BAM - Meta: $meta, BAM: $bam" }

        // SORTBAM required for perform_qcs, splitted_original_sample and split_by_chrom
        //if (params.perform_qcs | params.splitted_original_sample | params.split_by_chrom){
        SORTBAM(ALIGNRAWBAM.out.bam)
        bam_to_group = SORTBAM.out.bam.join(SORTBAM.out.csi)
        //}

        // truncate BAM to keep only the reads that are on target
        // TODO
        // see how BAMFILTERREADS requires the BAM file sorted....

        if (params.perform_qcs) {
            def qc_targets = params.targetsfile ?: []
            QUALIMAPQCRAW(SORTBAM.out.bam, ch_targetsfile)
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCRAW.out.results.map{it[1]}.collect())
        }

    
        if (params.splitted_original_sample){
            // Group BAMs by original sample name
            SORTBAM.out.bam
            .map { meta, bam -> 
                def sample = meta.sample
                tuple(sample, meta, bam)
            }
            .groupTuple(by: 0)
            .map { sample, metas, bams -> 
                def new_meta = metas[0].clone()
                new_meta.id = sample
                tuple(new_meta, bams)
            }
            .set { grouped_bams }

            // Merge BAMs for each sample
            MERGEBAM(grouped_bams)

            bam_to_group = MERGEBAM.out.bam_bai
        }

        if (params.split_by_chrom) {
            // Split merged BAMs by chromosome
            SPLITBAMCHROM(bam_to_group)

            def process_bams = { meta, bams ->
                bams.sort { it.name }.collect { bam ->
                    def new_meta = meta.clone()
                    // Extract chromosome from the last segment before .bam extension
                    // Expected format: <sample>.<chrom>.bam or <sample>_<chrom>.bam
                    // This prevents matching "chr" within the sample name itself
                    def chrom_match = bam.name =~ /\.(chr[^.]+)\.bam$/ 
                    def chrom = chrom_match ? chrom_match[0][1] : bam.name.replaceAll(/\.bam$/, '').replaceAll(/^.*[._]/, '')  
                    new_meta.id = "${meta.id}_${chrom}"  
                    [new_meta, bam]
                }
            }

            bam_to_group = SPLITBAMCHROM.out.chrom_bams
                .map { meta, bams -> process_bams(meta, bams) }
                .flatMap { it }
                .toSortedList { a, b -> a[0].id <=> b[0].id }
                .flatMap { it }
                .concat(
                    SPLITBAMCHROM.out.unknown_bam
                        .map { meta, bam -> 
                            def new_meta = meta.clone()
                            new_meta.id = "${meta.id}_unknown"
                            [new_meta, bam]
                        }
                )
                .map { meta, bam -> [meta, bam] } // Ensure correct structure
        }else {
            // The BAI index is dropped here because downstream processes only require the BAM file and its metadata.  
            def drop_bai_index = { meta, bam, bai -> tuple(meta, bam) }  
            bam_to_group = bam_to_group.map(drop_bai_index)  
        }

        //Final SORT for GROUPREADSBYUMI (TEMPLATE-COORDINATE SORTED)
        SORTBAMCLEAN(bam_to_group)
        non_duplex_bams = SORTBAMCLEAN.out.bam
    }
    //
    // Run fgbio Duplex consensus pipeline
    //

    if (params.step in ['mapping', 'groupreadsbyumi']) {

        // ASSIGN bam_to_group = to our input bam
        if (params.step == 'groupreadsbyumi') {
            non_duplex_bams = INPUT_CHECK.out.reads
        }

        // MODULE: Run fgbio GroupReadsByUmi
        // requires input template coordinate sorted
        GROUPREADSBYUMIDUPLEX(non_duplex_bams, "Paired")
        ch_multiqc_files = ch_multiqc_files.mix(GROUPREADSBYUMIDUPLEX.out.histogram.map{it[1]}.collect())

        // MODULE: Run fgbio CollecDuplexSeqMetrics
        COLLECTDUPLEXSEQMETRICS(GROUPREADSBYUMIDUPLEX.out.bam, [])
        
        // Extract family_sizes file directly from dedicated output
        family_sizes_metrics = COLLECTDUPLEXSEQMETRICS.out.family_sizes
        
        // When split_by_chrom is enabled, aggregate chromosome-specific files by sample
        if (params.split_by_chrom) {
            family_sizes_metrics = family_sizes_metrics
                .map { meta, file -> 
                    // Extract original sample name (remove chromosome suffix like "_chr1", "_chr2", "_unknown")
                    def original_sample = meta.sample ?: meta.id.replaceAll(/_(chr[^_]+|unknown)$/, '')
                    tuple(original_sample, meta, file)
                }
                .groupTuple(by: 0)  // Group by original sample name
                .map { sample, metas, files -> 
                    // Create new meta with original sample name
                    def new_meta = metas[0].clone()
                    new_meta.id = sample
                    new_meta.sample = sample
                    tuple(new_meta, files)
                }
        }
        
        // Plot the family size metrics
        FAMILYMETRICS(family_sizes_metrics)

        FAMILYMETRICS.out.sample_data.map{it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)
        FAMILYMETRICS.out.curve_data.map{it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetrics", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CollecDuplexSeqMetrics only on target
        COLLECTDUPLEXSEQMETRICSONTARGET(GROUPREADSBYUMIDUPLEX.out.bam, BEDTOINTERVAL.out.interval_list.first().map{it[1]} )
        
        // Extract family_sizes file directly from dedicated output
        family_sizes_metrics_ontarget = COLLECTDUPLEXSEQMETRICSONTARGET.out.family_sizes
        
        // When split_by_chrom is enabled, aggregate chromosome-specific files by sample
        if (params.split_by_chrom) {
            family_sizes_metrics_ontarget = family_sizes_metrics_ontarget
                .map { meta, file -> 
                    // Extract original sample name (remove chromosome suffix)
                    def original_sample = meta.sample ?: meta.id.replaceAll(/_(chr[^_]+|unknown)$/, '')
                    tuple(original_sample, meta, file)
                }
                .groupTuple(by: 0)  // Group by original sample name
                .map { sample, metas, files -> 
                    // Create new meta with original sample name
                    def new_meta = metas[0].clone()
                    new_meta.id = sample
                    new_meta.sample = sample
                    tuple(new_meta, files)
                }
        }

        // Plot the family size metrics
        FAMILYMETRICSONTARGET(family_sizes_metrics_ontarget)

        FAMILYMETRICSONTARGET.out.sample_data.map{it[1]}.collectFile(name: "metrics_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)
        FAMILYMETRICSONTARGET.out.curve_data.map{it[1]}.collectFile(name: "curves_summary.tsv", storeDir:"${params.outdir}/familymetricsontarget", skip: 1, keepHeader: true)


        // MODULE: Run fgbio CallDuplexConsensusReads
        CALLDUPLEXCONSENSUSREADS(GROUPREADSBYUMIDUPLEX.out.bam)
        
    }
    
    if (params.step in ['mapping', 'groupreadsbyumi', 'unmapped_consensus']) {

        if (params.step == 'unmapped_consensus') {
            UNMAPBAM(INPUT_CHECK.out.reads)
            called_consensus = UNMAPBAM.out.bam
        } else {
            called_consensus = CALLDUPLEXCONSENSUSREADS.out.bam
        }

        // MODULE: Align with bwa mem
        ALIGNDUPLEXCONSENSUSBAM(called_consensus, ch_ref_index_dir, false)

        SORTBAMALLMOLECULES(ALIGNDUPLEXCONSENSUSBAM.out.bam)

        if (params.split_by_chrom) {

            // Group BAMs by original sample name
            SORTBAMALLMOLECULES.out.bam
            .map { meta, bam -> 
                def sample = meta.sample
                tuple(sample, meta, bam)
            }
            .groupTuple(by: 0)
            .map { sample, metas, bams -> 
                def new_meta = metas[0].clone()
                new_meta.id = sample
                new_meta.sample = sample
                tuple(new_meta, bams)
            }
            .set { bam_n_index_all_molecules }

            // Merge BAMs for each sample
            MERGEBAMDUPLEX(bam_n_index_all_molecules)
            preASMINUSXSDUPLEXbams = MERGEBAMDUPLEX.out.bam_bai
        } else {
            SORTBAMALLMOLECULES.out.bam
                .join( SORTBAMALLMOLECULES.out.csi )
                .set { preASMINUSXSDUPLEXbams }
        }
    }

    if (params.step in ['mapping', 'groupreadsbyumi', 'unmapped_consensus', 'allmoleculesfile']) {

        if (params.step == 'allmoleculesfile') {
            preASMINUSXSDUPLEXbams = INPUT_CHECK.out.reads
        }

        ASMINUSXSDUPLEX(preASMINUSXSDUPLEXbams)
        SAMTOOLSFILTERALLMOLECULES(ASMINUSXSDUPLEX.out.bam)
        SORTBAMAMFILTERED(SAMTOOLSFILTERALLMOLECULES.out.bam)

        duplex_filtered_init_bam = SORTBAMAMFILTERED.out.bam

        ASMINUSXSDUPLEX.out.discarded_bam.map{[it[0], ch_targetsfile, it[1]]}.set { discarded_bam_targeted }
        DISCARDEDCOVERAGETARGETED(discarded_bam_targeted, [])

        ASMINUSXSDUPLEX.out.discarded_bam.map{[it[0], ch_global_exons_file, it[1]]}.set { discarded_bam }
        DISCARDEDCOVERAGEGLOBAL(discarded_bam, [])

    }

    if (params.step in ['mapping', 'groupreadsbyumi', 'unmapped_consensus', 'allmoleculesfile', 'filterconsensus']) {

        if (params.step == 'filterconsensus') {
            duplex_filtered_init_bam = INPUT_CHECK.out.reads
        }
   

        // Group by meta.parent_dna
        ch_grouped_bams = duplex_filtered_init_bam.map { meta, bam -> [['id' : meta.parent_dna], bam] }
                .groupTuple(by: 0)
                .filter { _meta, bams -> bams.size() >= 2 }
                .map { meta, bams -> [meta, bams.flatten()] }


        
        // Run the concatenation process
        MERGEBAMS(ch_grouped_bams)

        // Run sorting by query
        SORTBAMMERGED(MERGEBAMS.out.bam)

        duplex_filtered_bam = duplex_filtered_init_bam.mix(SORTBAMMERGED.out.bam)

        // store csv with all AM BAMs
        duplex_filtered_bam
            .map { meta, bam -> "sample,bam\n${meta.id},${params.outdir}/sortbamamfiltered/${bam.name}\n" }
            .collectFile(name: 'samplesheet_bam_filtered_inputs.csv', storeDir: "${params.outdir}/sortbamamfiltered", skip: 1, keepHeader: true)


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
            QUALIMAPQCDUPLEX(SORTBAMDUPLEXCONS.out.bam, ch_targetsfile)
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAPQCDUPLEX.out.results.map{it[1]}.collect())

            SORTBAMDUPLEXCONS.out.bam.map{it -> [it[0], ch_global_exons_file, it[1]]}.set { duplex_filt_bam_n_bed }
            COVERAGEGLOBAL(duplex_filt_bam_n_bed, [])

            // MULTIQCDUPLEX: Early report with accumulated QC files (no software versions yet)
            MULTIQCDUPLEX (
                ch_multiqc_files
                    .mix(Channel.from(ch_multiqc_config))
                    .mix(ch_multiqc_custom_config.collect().ifEmpty([]))
                    .collect()
            )
        }
    }

    if (params.step in ['mapping', 'groupreadsbyumi', 'unmapped_consensus', 'allmoleculesfile', 'filterconsensus', 'calling']) {
    
        // Initialize variables for calling step entry point
        if (params.step == 'calling') {
            cons_duplex_bam = INPUT_CHECK.out.reads
            bam_n_index_duplex_clean = INPUT_CHECK.out.reads
        }

        cons_duplex_bam.map{[it[0], it[1]]}
        .set{ cons_duplex_bam_only }

        // Compute depth of the consensus reads aligned to the genome
        COMPUTEDEPTH(cons_duplex_bam_only)

        CREATEBED(COMPUTEDEPTH.out.tsv)

        cons_duplex_bam
        .join( CREATEBED.out.bed )
        .set { cons_duplex_bam_bed }

        // Mutation calling for all reads using parallel VarDict subworkflow
        // Set vardict_chunks=1 for single-node processing if needed
        
        BAM_CALL_VARDICT_PARALLEL(
            cons_duplex_bam_bed,
            ch_ref_fasta,
            ch_ref_index_dir
        )

        // Postprocessing the BAM file to get exact coverage per position and allele
        //    also get the Ns per position
        RECOUNTMUTS(cons_duplex_bam,
                        bam_n_index_duplex_clean,
                        BAM_CALL_VARDICT_PARALLEL.out.vcf,
                        CREATEBED.out.bed,
                        ch_ref_fasta)

        if (params.annotate_mutations){
            VCFANNOTATE(BAM_CALL_VARDICT_PARALLEL.out.vcf,
                            ch_ref_fasta)
        }

        RECOUNTMUTS.out.somatic_vcf.map{it[1]}.set { mutation_files_duplex }
        SIGPROFPLOT(mutation_files_duplex.collect())

        RECOUNTMUTS.out.purvcf.map{it[1]}.set { mutation_files_pur_duplex }
        SIGPROFPLOTPUR(mutation_files_pur_duplex.collect())

        RECOUNTMUTS.out.pyrvcf.map{it[1]}.set { mutation_files_pyr_duplex }
        SIGPROFPLOTPYR(mutation_files_pyr_duplex.collect())

        // Generate deepCSA input example
        RECOUNTMUTS.out.filtered_vcf
        .join(cons_duplex_bam_only)
        .map { meta, vcf, bam -> "sample,vcf,bam\n${meta.id},${params.outdir}/mutations_vcf/${vcf.name},${params.outdir}/sortbamduplexcons/${bam.name}\n" }
        .collectFile(name: 'deepCSA_input_template.csv', storeDir: "${params.outdir}/pipeline_info", skip: 1, keepHeader: true)

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