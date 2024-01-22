
// Import required modules

include { BEDTOOLS_MERGE         as READJUSTREGIONS   } from '../../../modules/local/bedtools/merge/main'

include { SAMTOOLS_MPILEUP       as PILEUPBAM         } from '../../../modules/nf-core/samtools/mpileup/main'

include { NS_X_POSITION          as NSXPOSITION       } from '../../../modules/local/count_ns/main'

include { QUERY_TABIX            as QUERYTABIX        } from '../../../modules/local/tabix_mpileup/main'
include { PATCH_DEPTH            as PATCHDP           } from '../../../modules/local/patchdepth/main'

include { FILTER_LOW_COMPLEXITY  as FILTERLOWCOMPLEX  } from '../../../modules/local/filter/lowcomplexrep/main.nf'
include { FILTER_LOW_MAPPABILITY as FILTERLOWMAPPABLE } from '../../../modules/local/filter/lowmappability/main.nf'
include { FILTER_N_RICH          as FILTERNRICH       } from '../../../modules/local/filter/nrich/main.nf'

include { FILTERMUTATIONS      as FILTERVCF       } from '../../../modules/local/filtervcf/main'
include { MUTS_PER_POS         as MUTSPERPOS      } from '../../../modules/local/mutsperpos/main'



workflow RECOUNT_MUTS {

    take:

    bam_n_index              // channel: [mandatory] [ val(meta), path (bam), path (bamindex) ]
    vcf_file                 // channel: [mandatory] [ val(meta), path (vcf)]
    bed_file                 // channel: [mandatory] [ val(meta), path (intervals_file)]
    reference_fasta          // channel: [mandatory] path (reference_fasta)


    main:

    ch_versions = Channel.empty()

    vcf_file
    .join( bed_file )
    .set { ch_vcf_bed }

    READJUSTREGIONS(ch_vcf_bed)
    // These are the three main outputs
    // READJUSTREGIONS.out.vcf_bed
    // READJUSTREGIONS.out.vcf_bed_mut_ids
    // READJUSTREGIONS.out.regions_plus_variants_bed

    ch_versions = ch_versions.mix(READJUSTREGIONS.out.versions.first())

    // join the channel with the BAM file and the corresponding VCF
    // from the same samples together
    // SUPER IMPORTANT STEP
    bam_n_index
    .join( READJUSTREGIONS.out.regions_plus_variants_bed )
    .set { ch_bam_bai_bed }

    PILEUPBAM(ch_bam_bai_bed, reference_fasta)
    // This are the main outputs
    // PILEUPBAM.out.mpileup
    ch_versions = ch_versions.mix(PILEUPBAM.out.versions.first())


    NSXPOSITION(PILEUPBAM.out.mpileup)
    // This is the main output
    // NSXPOSITION.out.ns_per_pos
    ch_versions = ch_versions.mix(NSXPOSITION.out.versions.first())

    PILEUPBAM.out.mpileup
    .join( READJUSTREGIONS.out.vcf_bed )
    .set { ch_pileup_vcfbed }

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
    ch_versions = ch_versions.mix(PATCHDP.out.versions.first())
    
    // FINDMUTATED(ch_pileup_vcf)
    // FINDMUTATED.out.read_names
    // FINDMUTATED.out.tags
    // samtools view ../K_43_1_A_1_umi-grouped.bam -h -b -@ 9 -D MI:../K_43_1_A_1.tags > K_43_1_A_1.grouped.bam

    if (params.filter_mutations) {
        if (params.filter_human) {
            PATCHDP.out.patched_vcf
            .join( READJUSTREGIONS.out.vcf_bed_mut_ids )
            .set { ch_vcf_vcfbed }

            FILTERLOWCOMPLEX(ch_vcf_vcfbed)
            FILTERLOWMAPPABLE(FILTERLOWCOMPLEX.out.filtered_vcf_bed)
            
            FILTERLOWMAPPABLE.out.filtered_vcf
            .join(NSXPOSITION.out.ns_tsv)
            .set {ch_vcf_ns}
        } else {
            PATCHDP.out.patched_vcf
            .join(NSXPOSITION.out.ns_tsv)
            .set {ch_vcf_ns}
        }
        FILTERNRICH(ch_vcf_ns)
        output_vcf = FILTERNRICH.out.filtered_vcf
    } else {
        output_vcf = PATCHDP.out.patched_vcf
    }

    FILTERVCF(output_vcf)
    ch_versions = ch_versions.mix(FILTERVCF.out.versions.first())

    bam_n_index
    .join( FILTERVCF.out.vcf )
    .set { ch_bam_bai_vcf }
    
    MUTSPERPOS(ch_bam_bai_vcf)

    emit:

    ns_file        = NSXPOSITION.out.ns_tsv     // channel: [ val(meta), [ bed ], tbi ]
    versions       = ch_versions                // channel: [ versions.yml ]

    filtered_vcf   = output_vcf                 // channel: [ val(meta), [ vcf ] ]
    somatic_vcf    = FILTERVCF.out.vcf          // channel: [ val(meta), [ vcf ] ]


}
