//
// Run VEP to annotate VCF files
//

include { ENSEMBLVEP_VEP } from '../../../modules/nf-core/ensemblvep/vep/main'


workflow VCF_ANNOTATE_ENSEMBLVEP {
    take:
    vcf               // channel: [ val(meta), vcf ]
    fasta             //   value: fasta to use (optionnal)
    vep_genome        //   value: genome to use
    vep_species       //   value: species to use
    vep_cache_version //   value: cache version to use
    vep_cache         //    path: /path/to/vep/cache (optionnal)
    vep_extra_files   // channel: [ file1, file2...] (optionnal)

    main:
    

    ENSEMBLVEP_VEP(vcf, vep_genome, vep_species, vep_cache_version, vep_cache, fasta, vep_extra_files)


    emit:
    vcf      = ENSEMBLVEP_VEP.out.vcf      // channel: [ val(meta), vcf ]
    json     = ENSEMBLVEP_VEP.out.json     // channel: [ val(meta), json ]
    tab      = ENSEMBLVEP_VEP.out.tab      // channel: [ val(meta), tab ]
    reports  = ENSEMBLVEP_VEP.out.report   // channel: [ *.html ]
}