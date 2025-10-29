#!/usr/bin/env python3

import pysam
import click

def collect_retained_query_names(input_bam, threshold, tthreads):
    retained = set()
    with pysam.AlignmentFile(input_bam, "rb", threads=tthreads) as in_bam:
        for read in in_bam:
            if read.is_secondary:
                continue
            try:
                as_tag = read.get_tag("AS", with_value_type=False)
                xs_tag = read.get_tag("XS", with_value_type=False)
            except KeyError:
                continue
            if (as_tag - xs_tag) >= threshold:
                retained.add(read.query_name)
    return retained

def filter_bam(input_bam, output_bam, output_bam_discarded, retained_query_names, tthreads):
    with pysam.AlignmentFile(input_bam, "rb", threads=tthreads) as in_bam, \
         pysam.AlignmentFile(output_bam, "wb", header=in_bam.header, threads=tthreads) as out_bam, \
         pysam.AlignmentFile(output_bam_discarded, "wb", header=in_bam.header) as out_bam_discarded:
        for read in in_bam:
            if read.query_name in retained_query_names:
                out_bam.write(read)
            else:
                out_bam_discarded.write(read)

@click.command()
@click.option('--input_bam', type=click.Path(exists=True), required=True, help='Input BAM file (indexed).')
@click.option('--output_bam', type=click.Path(), required=True, help='Output BAM for retained reads.')
@click.option('--output_bam_discarded', type=click.Path(), required=True, help='Output BAM for discarded reads.')
@click.option('--threshold', type=int, required=True, help='AS-XS minimum difference threshold.')
@click.option('--tthreads', type=int, required=True, help='Number of threads to use when reading/writing BAMs.')
def main(input_bam, output_bam, output_bam_discarded, threshold, tthreads):
    """Filter BAM file based on AS minus XS threshold (efficient two-pass version)."""
    click.echo("First pass: collecting retained query names...")
    retained_query_names = collect_retained_query_names(input_bam, threshold, tthreads)
    click.echo(f"Collected {len(retained_query_names)} retained query names.")

    click.echo("Second pass: writing reads to output BAM files...")
    filter_bam(input_bam, output_bam, output_bam_discarded, retained_query_names, tthreads)
    click.echo("Done.")

if __name__ == "__main__":
    main()