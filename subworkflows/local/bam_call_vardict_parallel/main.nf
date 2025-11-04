/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW: Parallel VarDict variant calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Split BED file into chunks, call variants on each chunk in parallel jobs,
    then merge results back together.
----------------------------------------------------------------------------------------
*/

include { SPLIT_BED                 } from '../../../modules/local/calling_vardict/split_bed/main'
include { CALLING_VARDICT_CHUNK     } from '../../../modules/local/calling_vardict/calling_vardict_chunk/main'
include { MERGE_VARDICT_RESULTS     } from '../../../modules/local/calling_vardict/merge_vardict_results/main'

workflow BAM_CALL_VARDICT_PARALLEL {

    take:
    bam_bai_bed   // channel: [ val(meta), path(bam), path(bai), path(bed) ]
    fasta         // path: reference fasta
    fasta_dir     // path: directory containing reference fasta
    num_chunks    // val: number of chunks to split BED file into

    main:
    ch_versions = channel.empty()

    //
    // Step 1: Split BED file into chunks
    //
    SPLIT_BED(
        bam_bai_bed.map { meta, _bam, _bai, bed -> [meta, bed] },
        num_chunks
    )

    //
    // Step 2: Process each chunk in parallel
    // Flatten chunks and combine with BAM files for parallel processing
    //
    SPLIT_BED.out.chunks
        .transpose()  // Flatten list of chunks: [meta, [chunk1, chunk2, ...]] -> [meta, chunk1], [meta, chunk2], ...
        .map { meta, chunk -> 
            // Use meta.id as key for joining, keep original meta
            [meta.id, meta, chunk]
        }
        .combine(
            bam_bai_bed.map { meta, bam, bai, _bed -> [meta.id, bam, bai] },
            by: 0
        )
        .map { id, meta, chunk, bam, bai ->
            // Now add chunk name to meta for uniqueness
            def chunk_meta = meta.clone()
            chunk_meta.chunk_name = chunk.simpleName
            [chunk_meta, chunk, bam, bai]
        }
        .set { chunks_with_bam }

    CALLING_VARDICT_CHUNK(
        chunks_with_bam,
        fasta,
        fasta_dir
    )
    ch_versions = ch_versions.mix(CALLING_VARDICT_CHUNK.out.versions.first())

    //
    // Step 3: Merge all chunk results back together
    // Group chunks by original sample (without chunk_name)
    //
    CALLING_VARDICT_CHUNK.out.vcf
        .map { meta, vcf ->
            // Create clean meta without chunk_name for grouping
            def clean_meta = meta.clone()
            clean_meta.remove('chunk_name')
            [clean_meta, vcf]
        }
        .groupTuple()  // Group all chunks from same sample
        .set { grouped_vcf_chunks }

    // Also group raw TSV chunks by sample
    CALLING_VARDICT_CHUNK.out.raw
        .map { meta, raw ->
            def clean_meta = meta.clone()
            clean_meta.remove('chunk_name')
            [clean_meta, raw]
        }
        .groupTuple()
        .set { grouped_raw_chunks }

    // Join grouped VCF chunks and RAW chunks by meta
    grouped_vcf_chunks
        .join(grouped_raw_chunks)
        .set { grouped_vcf_and_raw }

    MERGE_VARDICT_RESULTS(grouped_vcf_and_raw)

    emit:
    vcf         = MERGE_VARDICT_RESULTS.out.vcf         // channel: [ val(meta), path(vcf) ]
    versions    = ch_versions                            // channel: [ path(versions.yml) ]
}
