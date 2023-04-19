
// Import required modules

include { FGBIO_FILTERBAM      as FGSELECTREADS   } from '../../../modules/local/fgbio/filterbam/main'

include { SAMTOOLS_SORT        as SORTBAMFIXED    } from '../../../modules/nf-core/samtools/sort/main'

include { SAMTOOLS_VIEW        as FILTERBAM       } from '../../../modules/nf-core/samtools/view/main'



workflow BAM_FILTER_READS {

    take:

    bam_n_index              // channel: [mandatory] [ val(meta), path (bam), path (bamindex) ]
    intervals_file           // channel: [mandatory] path (intervals_file)


    main:

    ch_versions = Channel.empty()

    FGSELECTREADS(bam_n_index, intervals_file)
    ch_versions = ch_versions.mix(FGSELECTREADS.out.versions.first())

    // join the created channel with the filtered reads to have
    // the ones from the same samples together
    // SUPER IMPORTANT STEP
    bam_n_index
    .join( FGSELECTREADS.out.read_names )
    .set { ch_bam_bai_reads }

    FILTERBAM(ch_bam_bai_reads,  [])
    ch_versions = ch_versions.mix(FILTERBAM.out.versions.first())

    SORTBAMFIXED(FILTERBAM.out.bam)
    ch_versions = ch_versions.mix(SORTBAMFIXED.out.versions.first())


    emit:

    bam        = SORTBAMFIXED.out.bam     // channel: [ val(meta), [ bam ] ]
    versions   = ch_versions              // channel: [ versions.yml ]
    
}
