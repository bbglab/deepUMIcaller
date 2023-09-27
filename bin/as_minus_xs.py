#!/usr/bin/env python

import pysam
import argparse

def filter_bam(input_bam, output_bam, threshold):
    with pysam.AlignmentFile(input_bam, "rb") as in_bam:
        header = in_bam.header
        with pysam.AlignmentFile(output_bam, "wb", header=header) as out_bam:
            for read in in_bam:
                if read.is_secondary:
                    # Skip secondary alignments
                    continue
                
                as_tag = read.get_tag("AS", with_value_type=False)
                xs_tag = read.get_tag("XS", with_value_type=False)
                
                if as_tag is not None and xs_tag is not None:
                    as_minus_xs = as_tag - xs_tag
                    if as_minus_xs >= threshold:
                        out_bam.write(read)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM file based on AS minus XS threshold.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("threshold", type=int, help="Threshold for AS minus XS")
    args = parser.parse_args()

    filter_bam(args.input_bam, args.output_bam, args.threshold)

if __name__ == "__main__":
    main()
