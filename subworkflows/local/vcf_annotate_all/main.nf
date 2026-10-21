//
// ANNOTATION
//

include { VCF_ANNOTATE_ENSEMBLVEP                       } from '../../nf-core/vcf_annotate_ensemblvep/main.nf'


workflow VCF_ANNOTATE_ALL {
    take:
    vcf          // channel: [ val(meta), vcf ]
    fasta


    main:
    reports = channel.empty()
    vcf_ann = channel.empty()
    tab_ann = channel.empty()
    json_ann = channel.empty()

    vep_extra_files = []


    VCF_ANNOTATE_ENSEMBLVEP(vcf, fasta,
                            params.vep_genome,
                            params.vep_species,
                            params.vep_cache_version,
                            params.vep_cache, vep_extra_files)


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