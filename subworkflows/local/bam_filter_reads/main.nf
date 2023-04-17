
// Import required modules

include { FGBIO_FILTERBAM      as FGSELECTREADS   } from '../../../modules/local/fgbio/filterbam/main'

include { SAMTOOLS_SORT        as SORTBAMFIXED    } from '../../../modules/nf-core/samtools/sort/main'

include { SAMTOOLS_VIEW        as FILTERBAM       } from '../../../modules/nf-core/samtools/view/main'



workflow BAM_FILTER_READS {

    take:

    bam                      // channel: [mandatory] [ val(meta), path (bam) ]
    bamindex                 // channel: [mandatory] path (bamindex)
    intervals_file           // channel: [mandatory] path (intervals_file)


    main:

    ch_versions = Channel.empty()

    FGSELECTREADS(bam, intervals_file)
    ch_versions = ch_versions.mix(FGSELECTREADS.out.versions.first())

    FILTERBAM(bam, bamindex,  [], FGSELECTREADS.out.read_names.map{it -> it[1]} )
    ch_versions = ch_versions.mix(FILTERBAM.out.versions.first())

    SORTBAMFIXED(FILTERBAM.out.bam)
    ch_versions = ch_versions.mix(SORTBAMFIXED.out.versions.first())


    emit:

    bam        = SORTBAMFIXED.out.bam     // channel: [ val(meta), [ bam ] ]
    versions   = ch_versions              // channel: [ versions.yml ]
    
}
