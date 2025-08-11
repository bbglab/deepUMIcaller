
// Import required modules

include { BEDTOOLS_MERGE         as READJUSTREGIONS   } from '../../../modules/local/bedtools/merge/main'

include { SAMTOOLS_MPILEUP       as PILEUPBAM         } from '../../../modules/nf-core/samtools/mpileup/main'
include { SAMTOOLS_MPILEUP       as PILEUPBAMALL      } from '../../../modules/nf-core/samtools/mpileup/main'

include { NS_X_POSITION          as NSXPOSITION       } from '../../../modules/local/count_ns/main'

include { QUERY_TABIX            as QUERYTABIX        } from '../../../modules/local/tabix_mpileup/main'
include { PATCH_DEPTH            as PATCHDP           } from '../../../modules/local/patchdepth/main'
include { PATCH_DEPTH            as PATCHDPALL        } from '../../../modules/local/patchdepth/main'

include { FILTER_LOW_COMPLEXITY  as FILTERLOWCOMPLEX  } from '../../../modules/local/filter/lowcomplexrep/main.nf'
include { FILTER_LOW_MAPPABILITY as FILTERLOWMAPPABLE } from '../../../modules/local/filter/lowmappability/main.nf'
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

    low_complex_filter = params.low_complex_file ? Channel.fromPath( params.low_complex_file, checkIfExists: true).first() : Channel.fromPath(params.input)
    low_mappability_filter = params.low_mappability_file ? Channel.fromPath( params.low_mappability_file, checkIfExists: true).first() : Channel.fromPath(params.input)

    

    vcf_file
    .join( bed_file )
    .set { ch_vcf_bed }

    READJUSTREGIONS(ch_vcf_bed)

    

    // join the channel with the BAM file and the corresponding VCF
    // from the same samples together
    bam_n_index
    .join( READJUSTREGIONS.out.regions_plus_variants_bed )
    .set { ch_bam_bai_bed }

    PILEUPBAM(ch_bam_bai_bed, reference_fasta)
    


    NSXPOSITION(PILEUPBAM.out.mpileup)
    

    PILEUPBAM.out.mpileup
    .join( READJUSTREGIONS.out.vcf_bed )
    .set { ch_pileup_vcfbed }


    // join the channel with the BAM file and the corresponding VCF
    // from the same samples together
    bam_n_index_all_mol
    .join( READJUSTREGIONS.out.vcf_bed )
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


    if (params.filter_regions) {
        PATCHDPALL.out.patched_vcf
        .join( READJUSTREGIONS.out.vcf_bed_mut_ids )
        .set { ch_vcf_vcfbed }

        FILTERLOWCOMPLEX(ch_vcf_vcfbed, low_complex_filter)
        FILTERLOWMAPPABLE(FILTERLOWCOMPLEX.out.filtered_vcf_bed, low_mappability_filter)
        
        FILTERLOWMAPPABLE.out.filtered_vcf
        .join(NSXPOSITION.out.ns_tsv)
        .set {ch_vcf_ns}
    } else {
        PATCHDPALL.out.patched_vcf
        .join(NSXPOSITION.out.ns_tsv)
        .set {ch_vcf_ns}
    }
    FILTERNRICH(ch_vcf_ns)
    output_vcf = FILTERNRICH.out.filtered_vcf

    FILTERVCFSOMATIC(output_vcf)
    

    FILTERVCFPLOT(output_vcf)
    

    bam_n_index
    .join( FILTERVCFPLOT.out.vcf )
    .set { ch_bam_bai_vcf }
    
    MUTSPERPOS(ch_bam_bai_vcf)
    MUTSPERPOS.out.positions_csv.map{ it[1] }.collect().set{ mutations_position_csv }

    COHORTMUTSPERPOS(mutations_position_csv)

    emit:

    ns_file         = NSXPOSITION.out.ns_tsv     // channel: [ val(meta), [ bed ], tbi ]

    filtered_vcf    = output_vcf                 // channel: [ val(meta), [ vcf ] ]
    somatic_vcf     = FILTERVCFSOMATIC.out.vcf   // channel: [ val(meta), [ vcf ] ]

    purvcf          = FILTERVCFSOMATIC.out.pur_vcf
    pyrvcf          = FILTERVCFSOMATIC.out.pyr_vcf

}
