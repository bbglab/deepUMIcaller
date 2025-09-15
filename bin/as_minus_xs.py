#!/usr/bin/env python

import pysam
import argparse

def filter_bam(input_bam, output_bam, output_bam_discarded, threshold, tthreads):
    with pysam.AlignmentFile(input_bam, "rb", threads=tthreads) as in_bam:
        with pysam.AlignmentFile(output_bam, "wb", header=in_bam.header, threads=tthreads) as out_bam:
            with pysam.AlignmentFile(output_bam_discarded, "wb", header=in_bam.header) as out_bam_discarded:
                # Create a set to store query names
                query_names = set()
                # Dictionary to store discarded reads
                reads_discarded = {}

                for read in in_bam:
                    if read.is_secondary:
                        # Skip secondary alignments
                        continue
                    
                    try:
                        as_tag = read.get_tag("AS", with_value_type=False)
                        xs_tag = read.get_tag("XS", with_value_type=False)
                    except KeyError:
                        continue  # Skip reads without AS/XS tags
                    
                    # Check if the query name is not in the set of (already checked) query names
                    if read.query_name not in query_names:
                        if as_tag is not None and xs_tag is not None:
                            as_minus_xs = as_tag - xs_tag
                            if as_minus_xs >= threshold:
                                # Add the query name to the set of query names
                                query_names.add(read.query_name)
                                out_bam.write(read)
                                continue

                            # Store potentially discarded reads in a dictionary
                            elif read.query_name not in reads_discarded.keys():
                                reads_discarded[read.query_name] = set()

                            reads_discarded[read.query_name].add(read)
                        
                    else:
                        # Write retained reads to the output BAM file
                        out_bam.write(read)

                # Write potentially discarded reads to the appropriate output BAM file
                for key, value in reads_discarded.items():
                    if key in query_names:
                        # If the query name is retained, write the discarded reads to the output BAM file for retained reads
                        for read in value:
                            out_bam.write(read)
                    else:
                        # If the query name is discarded, write the discarded reads to the output BAM file for discarded reads
                        for read in value:
                            out_bam_discarded.write(read)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM file based on AS minus XS threshold.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file")
    parser.add_argument("output_bam_discarded", help="Output discarded BAM file")
    parser.add_argument("threshold", type=int, help="Threshold for AS minus XS")
    parser.add_argument("tthreads", type=int, help="threads")
    args = parser.parse_args()

    filter_bam(args.input_bam, args.output_bam, args.output_bam_discarded, args.threshold, args.tthreads)

if __name__ == "__main__":
    main()
