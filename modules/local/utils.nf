
def process_bams(meta, bam) {

    def new_meta = meta.clone()
    def filename = bam.name

    // Extract chromosome from the last segment before .bam extension
    // Expected format: <sample>.<chrom>.bam or <sample>_<chrom>.bam
    // This prevents matching "chr" within the sample name itself
    def chrom_match = filename =~ /\.(chr[^.]+)\.bam$/ 
    def chrom = chrom_match ? chrom_match[0][1] : filename.replaceAll(/\.bam$/, '').replaceAll(/^.*[._]/, '')  
    new_meta.id = "${meta.id}_${chrom}"

    return [new_meta, bam]
}