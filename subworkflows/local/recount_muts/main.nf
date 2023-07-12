
// Import required modules

include { BEDTOOLS_MERGE       as READJUSTREGIONS } from '../../../modules/local/bedtools/merge/main'

include { SAMTOOLS_MPILEUP     as PILEUPBAM       } from '../../../modules/nf-core/samtools/mpileup/main'

include { PATCH_DEPTH          as PATCHDP         } from '../../../modules/local/patchdepth/main'

include { NS_X_POSITION        as NSXPOSITION     } from '../../../modules/local/count_ns/main'


workflow RECOUNT_MUTS {

    take:

    bam_n_index              // channel: [mandatory] [ val(meta), path (bam), path (bamindex) ]
    vcf_file                 // channel: [mandatory] [ val(meta), path (vcf)]
    bed_file                 // channel: [mandatory] path (intervals_file)
    reference_fasta          // channel: [mandatory] path (reference_fasta)


    main:

    ch_versions = Channel.empty()

    READJUSTREGIONS(vcf_file, bed_file)
    // These are the two main outputs
    // READJUSTREGIONS.out.vcf_bed
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


    NSXPOSITION(PILEUPBAM.out.mpileup) // we could consider tabixing this file
    // This is the main output
    // NSXPOSITION.out.ns_per_pos
    ch_versions = ch_versions.mix(NSXPOSITION.out.versions.first())

    PILEUPBAM.out.mpileup
    .join( READJUSTREGIONS.out.vcf_bed )
    .set { ch_pileup_vcfbed }

    ch_pileup_vcfbed
    .join( vcf_file )
    .set { ch_pileup_vcfbed_vcf }


    PATCHDP(ch_pileup_vcfbed_vcf)
    // This is the main output
    // PATCHDP.out.patched_vcf
    // think well which is the best way to output this information, if a VCF or a TSV with only the updated depths or what.
    // also think whether it makes sense to remove strand bias flags from the VCF file
    //   maybe it makes 
    ch_versions = ch_versions.mix(PATCHDP.out.versions.first())


    emit:

    ns_file        = NSXPOSITION.out.ns_tsv     // channel: [ val(meta), [ bed ], tbi ]
    corrected_vcf  = PATCHDP.out.patched_vcf    // channel: [ val(meta), [ bam ] ]
    versions       = ch_versions                // channel: [ versions.yml ]

}
