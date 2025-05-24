process SAMTOOLS_FILTER {
    tag "$meta.id"
    label 'process_medium'
    label 'process_medium_high_memory'

    container 'quay.io/pacbio/samtools:1.17'
    // conda "bioconda::samtools=1.16.1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
    //     'biocontainers/samtools:1.16.1--h6899075_1' }"

    // conda "bioconda::samtools=1.16.1"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1' :
    //     'biocontainers/samtools:1.16.1--h6899075_1' }"
    // container 'mgibio/samtools-cwl'

    input:
    tuple val(meta), path(bam)
    val threshold

    output:
    tuple val(meta), path("*.0x2.AS-XS.bam") , emit: bam ,    optional: true
    tuple val(meta), path("*.cram"), emit: cram,    optional: true
    tuple val(meta), path("*.sam") , emit: sam ,    optional: true
    tuple val(meta), path("*.bai") , emit: bai ,    optional: true
    tuple val(meta), path("*.csi") , emit: csi ,    optional: true
    tuple val(meta), path("*.crai"), emit: crai,    optional: true
    path  "versions.yml"           , emit: versions


    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = args.contains("--output-fmt sam") ? "sam" :
                    args.contains("--output-fmt bam") ? "bam" :
                    args.contains("--output-fmt cram") ? "cram" :
                    bam.getExtension()
    if ("$bam" == "${prefix}.${file_type}") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools view ${bam} -b -h -f 0x2 > ${prefix}.filtered.0x2.bam

    cat <<EOF > filter_bam.cpp
    #include <iostream>
    #include <cstring>
    #include <htslib/sam.h>

    int main(int argc, char *argv[]) {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " <input.bam> <output.bam>" << std::endl;
            return 1;
        }

        const char *input_bam = argv[1];
        const char *output_bam = argv[2];
        samFile *in = sam_open(input_bam, "r");

        if (in == NULL) {
            std::cerr << "Error opening input BAM file: " << input_bam << std::endl;
            return 1;
        }

        bam_hdr_t *header = sam_hdr_read(in);
        samFile *out = sam_open(output_bam, "wb");

        if (out == NULL) {
            std::cerr << "Error opening output BAM file: " << output_bam << std::endl;
            return 1;
        }

        sam_hdr_write(out, header) __attribute__((warn_unused_result)); // Suppress warning
        bam1_t *record = bam_init1();

        while (sam_read1(in, header, record) >= 0) {
            uint32_t as = 0, xs = 0;

            // Extract AS and XS fields from the BAM record
            uint8_t *as_ptr = bam_aux_get(record, "AS");
            uint8_t *xs_ptr = bam_aux_get(record, "XS");

            if (as_ptr != NULL) as = bam_aux2i(as_ptr);
            if (xs_ptr != NULL) xs = bam_aux2i(xs_ptr);

            int as_minus_xs = as - xs;

            if (as_minus_xs >= ${threshold}) {
                sam_write1(out, header, record) __attribute__((warn_unused_result)); // Suppress warning
            }
        }

        bam_destroy1(record);
        bam_hdr_destroy(header);
        sam_close(in);
        sam_close(out);

        return 0;
    }
    EOF

    g++ -o filter_bam filter_bam.cpp -Wl,-rpath,/usr/local/lib -I /usr/local/include/ -L /usr/local/lib -lhts;

    ./filter_bam ${prefix}.filtered.0x2.bam ${prefix}.filtered.0x2.AS-XS.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam
    touch ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
    // samtools \\
    //     view \\
    //     --threads ${task.cpus-1} \\
    //     $args \\
    //     ${reference} \\
    //     ${readnames} \\
    //     -o ${prefix}.${file_type} \\
    //     $input \\
    //     $args2
