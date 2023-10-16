#!/usr/bin/env python

import pysam
import argparse

def filter_bam(input_bam, output_bam, threshold):
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        header = in_bam.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            written_reads = set()  # Set to store read IDs that have been written
            for read in in_bam:
                # Check if the read ID is in the set of written reads
                if read.query_name in written_reads:
                    continue
                
                # Check the AS and XS tags
                as_tag = read.get_tag("AS")
                xs_tag = read.get_tag("XS")
                
                # Check the mapping quality of both the read and its pair
                if (
                    as_tag is not None
                    and xs_tag is not None
                    and read.mapping_quality >= threshold
                ):
                    as_minus_xs = as_tag - xs_tag
                    if as_minus_xs >= threshold:
                        # Write the read and its pair to the output BAM file
                        out_bam.write(read)
                        written_reads.add(read.query_name)
                        if read.is_paired and read.is_proper_pair:
                            pair = in_bam.mate(read)
                            pair_as_tag = pair.get_tag("AS")
                            pair_xs_tag = pair.get_tag("XS")
                            if (
                                pair_as_tag is not None
                                and pair_xs_tag is not None
                            ):
                                pair_as_minus_xs = pair_as_tag - pair_xs_tag
                                if pair_as_minus_xs >= threshold:
                                    out_bam.write(pair)
                                    written_reads.add(pair.query_name)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM file based on AS minus XS threshold and minimum quality.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("threshold", type=int, help="Threshold for AS minus XS and minimum quality")
    args = parser.parse_args()

    filter_bam(args.input_bam, args.output_bam, args.threshold)

if __name__ == "__main__":
    main()
