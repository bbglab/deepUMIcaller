//
// ANNOTATION
//

include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main.nf'


workflow VCF_ANNOTATE_ALL {
    take:
    vcf          // channel: [ val(meta), vcf ]
    fasta
    vep_genome
    vep_species
    vep_cache_version
    vep_cache
    vep_extra_files

    main:
    reports = Channel.empty()
    vcf_ann = Channel.empty()
    tab_ann = Channel.empty()
    json_ann = Channel.empty()

    // if (tools.split(',').contains('merge') || tools.split(',').contains('snpeff')) {
    //     VCF_ANNOTATE_SNPEFF(vcf, snpeff_db, snpeff_cache)

    //     reports = reports.mix(VCF_ANNOTATE_SNPEFF.out.reports)
    //     vcf_ann = vcf_ann.mix(VCF_ANNOTATE_SNPEFF.out.vcf_tbi)
    //     versions = versions.mix(VCF_ANNOTATE_SNPEFF.out.versions)
    // }

    // if (tools.split(',').contains('merge')) {
    //     vcf_ann_for_merge = VCF_ANNOTATE_SNPEFF.out.vcf_tbi.map{ meta, vcf, tbi -> [ meta, vcf ] }
    //     VCF_ANNOTATE_MERGE(vcf_ann_for_merge, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

    //     reports = reports.mix(VCF_ANNOTATE_MERGE.out.reports)
    //     vcf_ann = vcf_ann.mix(VCF_ANNOTATE_MERGE.out.vcf_tbi)
    //     versions = versions.mix(VCF_ANNOTATE_MERGE.out.versions)
    // }

    VCF_ANNOTATE_ENSEMBLVEP(vcf, fasta, vep_genome, vep_species, vep_cache_version, vep_cache, vep_extra_files)

    reports = reports.mix(VCF_ANNOTATE_ENSEMBLVEP.out.reports)
    vcf_ann = vcf_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.vcf)
    tab_ann = tab_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.tab)
    json_ann = json_ann.mix(VCF_ANNOTATE_ENSEMBLVEP.out.json)

    emit:
    vcf_ann      // channel: [ val(meta), vcf.gz, vcf.gz.tbi ]
    tab_ann
    json_ann
    reports      //    path: *.html
}