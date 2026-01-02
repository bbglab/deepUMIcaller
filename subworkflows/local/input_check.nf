//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    step

    main:
    SAMPLESHEET_CHECK ( samplesheet, step )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { it -> create_input_channel(it,step) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row, step) {
    
    // create meta map
    def meta = [:]

    meta.id             = row.id
    meta.sample         = row.sample

    // Set parent_dna to row.parent_dna if present and not empty, else use id to avoid unintended merging
    meta.parent_dna     = (row.containsKey('parent_dna') && row.parent_dna) ? row.parent_dna : row.id

    // add path(s) of the fastq file(s) to the meta map
    def input_meta = []
    if (step=='mapping'){
        meta.read_structure = row.read_structure
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        input_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }
    // For all the calling, only BAMs/indexes are required
    else if (step=='calling' || step=='allmoleculesfile'){
        if (!file(row.duplexbam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Duplex BAM file does not exist!\n${row.duplexbam}"
        }
        if (!file(row.csi).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Csi file does not exist!\n${row.csi}"
        }
        input_meta = [ meta, [ file(row.duplexbam) ], [ file(row.csi) ] ]
    }

    // Only the bam is required to the meta map for the other steps (e.g. filterconsensus)
    else{
        if (!file(row.bam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> BAM file does not exist!\n${row.bam}"
        }
        input_meta = [ meta, [ file(row.bam) ]]
    }
    return input_meta

}
