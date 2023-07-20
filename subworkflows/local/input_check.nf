//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    step

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_input_channel(it,step) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_input_channel(LinkedHashMap row, step) {
    
    // create meta map
    def meta = [:]
    meta.id             = row.sample
    meta.read_structure = row.read_structure

    // add path(s) of the fastq file(s) to the meta map
    def input_meta = []
    if (step=='mapping'){
        if (!file(row.fastq_1).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
        }
        if (!file(row.fastq_2).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        input_meta = [ meta, [ file(row.fastq_1), file(row.fastq_2) ] ]
    }

    // For all the other steps, only bams/indexes are required
    else{
        if (!file(row.bam).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Bam file does not exist!\n${row.bam}"
        }
        if (!file(row.bai).exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Bai file does not exist!\n${row.bai}"
        }
        input_meta = [ meta, [ file(row.bam), file(row.bai) ] ]
    }
    return input_meta

}
