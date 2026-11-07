def process_bams(meta, bams) {
    def results = []
    bams.each { bam ->
        def new_meta = meta.clone()
        def filename = bam.name

        // Extract chromosome from the last segment before .bam extension
        // Expected format: <sample>.<chrom>.bam or <sample>_<chrom>.bam
        // This prevents matching "chr" within the sample name itself
        def chrom_match = filename =~ /\.(chr[^.]+)\.bam$/
        def chrom = chrom_match ? chrom_match[0][1] : filename.replaceAll(/\.bam$/, '').replaceAll(/^.*[._]/, '')
        new_meta.id = "${meta.id}_${chrom}"

        results << [new_meta, bam]
    }
    return results
}

def clean_chr_names(meta_file_pairs) {
    def sample_grouped_file_pairs
    sample_grouped_file_pairs = meta_file_pairs.map { meta, file -> 
                    // Extract original sample name (remove chromosome suffix)
                    def original_sample = meta.sample ?: meta.id.replaceAll(/_(chr[^_]+|unknown)$/, '')
                    tuple(original_sample, file)
                }
                .groupTuple(by: 0)  // Group by original sample name
                .map { sample, files -> 
                    // Create new meta with original sample name
                    def new_meta = [id: sample, sample: sample]
                    tuple(new_meta, files.sort { it -> it.name })
                }
    return sample_grouped_file_pairs
}