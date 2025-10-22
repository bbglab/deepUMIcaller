
// Import required modules

include { BEDTOOLS_MERGE         as READJUSTREGIONS_AMP   } from '../../../modules/local/bedtools/merge/main'
include { BEDTOOLS_MERGE         as READJUSTREGIONS_NOAMP   } from '../../../modules/local/bedtools/merge/main'

include { SAMTOOLS_MPILEUP       as PILEUPBAM         } from '../../../modules/nf-core/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP       as PILEUPBAMALL      } from '../../../modules/nf-core/samtools/mpileup/main'

include { NS_X_POSITION          as NSXPOSITION       } from '../../../modules/local/count_ns/main'

include { QUERY_TABIX            as QUERYTABIX        } from '../../../modules/local/tabix_mpileup/main'
include { PATCH_DEPTH            as PATCHDP           } from '../../../modules/local/patchdepth/main'
include { PATCH_DEPTH            as PATCHDPALL        } from '../../../modules/local/patchdepth/main'

include { FILTER_FROM_BED        as FILTERLOWMAPPABLE } from '../../../modules/local/filter/from_bed/main.nf'
include { FILTER_FROM_BED        as FILTERLOWCOMPLEX  } from '../../../modules/local/filter/from_bed/main.nf'
include { FILTER_FROM_BED        as FILTERNANOSEQSNP  } from '../../../modules/local/filter/from_bed/main.nf'
include { FILTER_FROM_BED        as FILTERNANOSEQNOISE} from '../../../modules/local/filter/from_bed/main.nf'
include { FILTER_N_RICH          as FILTERNRICH       } from '../../../modules/local/filter/nrich/main.nf'

include { FILTERMUTATIONS        as FILTERVCFSOMATIC  } from '../../../modules/local/filtervcf/main'
include { FILTERMUTATIONS        as FILTERVCFPLOT     } from '../../../modules/local/filtervcf/main'

include { MUTS_PER_POS           as MUTSPERPOS        } from '../../../modules/local/mutsperpos/compute/main'
include { SUMMARIZE_MUTS_PER_POS as COHORTMUTSPERPOS  } from '../../../modules/local/mutsperpos/summarize/main'



workflow RECOUNT_MUTS {

    take:

    bam_n_index              // channel: [mandatory] [ val(meta), path (bam), path (bamindex) ]
    bam_n_index_all_mol      // channel: [mandatory] [ val(meta), path (bam), path (bamindex) ]
    vcf_file                 // channel: [mandatory] [ val(meta), path (vcf)]
    bed_file                 // channel: [mandatory] [ val(meta), path (intervals_file)]
    reference_fasta          // channel: [mandatory] path (reference_fasta)


    main:

    low_complex_filter = params.low_complex_file ? Channel.fromPath(params.low_complex_file, checkIfExists: true).first() : Channel.empty()
    low_mappability_filter = params.low_mappability_file ? Channel.fromPath( params.low_mappability_file, checkIfExists: true).first() : Channel.empty()
    nanoseq_snp_filter = params.nanoseq_snp_file ? Channel.fromPath( params.nanoseq_snp_file, checkIfExists: true).first() : Channel.empty()
    nanoseq_noise_filter = params.nanoseq_noise_file ? Channel.fromPath( params.nanoseq_noise_file, checkIfExists: true).first() : Channel.empty()

    vcf_file
    .join( bed_file )
    .set { ch_vcf_bed }

    // With amplification
    readjust_amp_out = READJUSTREGIONS_AMP(ch_vcf_bed)
    ch_readjustregions_amp = readjust_amp_out.vcf_bed_mut_ids
    // Without amplification
    readjust_noamp_out = READJUSTREGIONS_NOAMP(ch_vcf_bed)
    ch_readjustregions_noamp = readjust_noamp_out.vcf_bed_mut_ids


    // join the channel with the BAM file and the corresponding VCF
    // from the same samples together
    bam_n_index
    .join( readjust_amp_out.regions_plus_variants_bed )
    .set { ch_bam_bai_bed }

    PILEUPBAM(ch_bam_bai_bed, reference_fasta)
    

    ch_nsxposition =  NSXPOSITION(PILEUPBAM.out.mpileup)
    

    PILEUPBAM.out.mpileup
    .join( readjust_amp_out.vcf_bed )
    .set { ch_pileup_vcfbed }


    // join the channel with the BAM file and the corresponding VCF
    // from the same samples together
    bam_n_index_all_mol
    .join( readjust_amp_out.vcf_bed )
    .set { ch_bamall_bai_bed }

    PILEUPBAMALL(ch_bamall_bai_bed, reference_fasta)
    
    QUERYTABIX(ch_pileup_vcfbed)

    QUERYTABIX.out.mutated_tsv
    .join( vcf_file )
    .set { ch_pileup_vcf }

    PATCHDP(ch_pileup_vcf)
    // This is the main output
    // PATCHDP.out.patched_vcf
    // think well which is the best way to output this information, if a VCF or a TSV with only the updated depths or what.
    // also think whether it makes sense to remove strand bias flags from the VCF file
    //   maybe it makes 
    

    PILEUPBAMALL.out.mpileup.map{ it -> [it[0], it[1]] }
    .join( PATCHDP.out.patched_vcf )
    .set{ ch_pileup_vcfpatched1 }

    PATCHDPALL(ch_pileup_vcfpatched1)


    def do_filter_regions = params.filter_regions
    def has_low_mappability_file = params.low_mappability_file
    def has_low_complex_file = params.low_complex_file

    // Validate if specie is human to use nanoseq filters
    def has_nanoseq_snp_file = (params.vep_species == "homo_sapiens") ? params.nanoseq_snp_file : null
    def has_nanoseq_noise_file = (params.vep_species == "homo_sapiens") ? params.nanoseq_noise_file : null

    if (params.vep_species != "homo_sapiens" && (params.nanoseq_snp_file || params.nanoseq_noise_file)) {
        log.warn "Nanoseq filters unset for other species than homo sapiens"
    }

    // First, always create ch_vcf_vcfbed
    PATCHDPALL.out.patched_vcf
        .join(ch_readjustregions_amp)
        .set { ch_vcf_vcfbed }
    
    def ch_vcf_current = ch_vcf_vcfbed

    // Apply low mappability filter if enabled
    ch_vcf_lowmappable = do_filter_regions && has_low_mappability_file ?
        FILTERLOWMAPPABLE(
            ch_vcf_current.combine(low_mappability_filter)
                .map { meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path -> 
                    [meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path, "low_mappability"]
                }
        ).filtered_vcf_bed :
        ch_vcf_current
    
    // Apply low complex filter if enabled
    ch_vcf_lowcomplex = do_filter_regions && has_low_complex_file ?
        FILTERLOWCOMPLEX(
            ch_vcf_lowmappable.combine(low_complex_filter)
                .map { meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path ->
                    [meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path, "low_complex"]
                }
        ).filtered_vcf_bed :
        ch_vcf_lowmappable
    
    // Join with readjusted regions
    ch_vcf_current_noamp = ch_vcf_lowcomplex
        .join(ch_readjustregions_noamp)
        .map { meta, vcf_file_path, vcf_derived_bed_path, vcf_derived_bed_noamp -> [meta, vcf_file_path, vcf_derived_bed_noamp] }

    // Apply nanoseq SNP and noise filters if enabled
    if (do_filter_regions && has_nanoseq_snp_file && has_nanoseq_noise_file) {
        ch_vcf_nanoseqsnp = FILTERNANOSEQSNP(
            ch_vcf_current_noamp.combine(nanoseq_snp_filter)
            .map { meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path ->
                    [meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path, "nanoseq_snp"]
            }
            ).filtered_vcf_bed
        ch_vcf_nanoseqnoise = FILTERNANOSEQNOISE(ch_vcf_nanoseqsnp.combine(nanoseq_noise_filter)
            .map { meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path ->
                [meta, vcf_file_path, vcf_derived_bed_path, mask_bed_path, "nanoseq_noise"]
            }).filtered_vcf_bed
        ch_vcf_final = ch_vcf_nanoseqnoise.join(ch_nsxposition.ns_tsv)
    } else {
        log.warn "No bed files provided for nanoseq filters, filter_regions=false or species is not homo sapiens; skipping bed-based region filters."
        ch_vcf_final = ch_vcf_lowcomplex.join(ch_nsxposition.ns_tsv)
    }

    // FILTER N RICH
    ch_vcf_final
    .map { tuple -> 
        def (meta, vcf_file_path, vcf_derived_bed_path, ns_position_file, ns_position_index) = tuple
        [meta, vcf_file_path, ns_position_file, ns_position_index]
    }
    .set { ch_vcf_nrich_input }

    FILTERNRICH(ch_vcf_nrich_input)
    output_vcf = FILTERNRICH.out.filtered_vcf

    // FILTER SOMATIC VARIANTS
    FILTERVCFSOMATIC(output_vcf)
    

    FILTERVCFPLOT(output_vcf)
    

    bam_n_index
    .join( FILTERVCFPLOT.out.vcf )
    .set { ch_bam_bai_vcf }
    
    MUTSPERPOS(ch_bam_bai_vcf)
    MUTSPERPOS.out.positions_csv.map{ it[1] }.collect().set{ mutations_position_csv }

    COHORTMUTSPERPOS(mutations_position_csv)

    emit:

    ns_file         = ch_nsxposition.ns_tsv     // channel: [ val(meta), [ bed ], tbi ]

    filtered_vcf    = output_vcf                 // channel: [ val(meta), [ vcf ] ]
    somatic_vcf     = FILTERVCFSOMATIC.out.vcf   // channel: [ val(meta), [ vcf ] ]

    purvcf          = FILTERVCFSOMATIC.out.pur_vcf
    pyrvcf          = FILTERVCFSOMATIC.out.pyr_vcf

}
