# BAM_CALL_VARDICT_PARALLEL

## Description

A subworkflow for parallel VarDict variant calling across multiple compute nodes. This subworkflow splits genomic regions into chunks, processes each chunk independently on separate nodes, and merges the results.

## Usage

```nextflow
include { BAM_CALL_VARDICT_PARALLEL } from './subworkflows/local/bam_call_vardict_parallel/main'

workflow {
    BAM_CALL_VARDICT_PARALLEL(
        bam_bai_bed_ch,    // channel: [ val(meta), path(bam), path(bai), path(bed) ]
        fasta,             // path: reference.fa
        fasta_dir          // path: directory containing reference.fa
    )
}
```

## Input

- `bam_bai_bed`: Channel containing:
  - `meta`: Map with sample metadata (e.g., `[id: 'sample1']`)
  - `bam`: Sorted BAM file
  - `bai`: BAM index file
  - `bed`: BED file with target regions
- `fasta`: Reference genome FASTA file
- `fasta_dir`: Directory containing the reference genome

**Note**: The number of chunks is controlled by the `params.vardict_chunks` parameter in `nextflow.config` (defaults to 2, or falls back to `params.max_cpus` if not set)

## Output

- `vcf`: Channel with filtered VCF files `[ val(meta), path(*.vcf) ]`
- `tsv`: Channel with merged, gzipped raw TSV files `[ val(meta), path(*.tsv.gz) ]`
- `genome_vcf`: Channel with compressed unfiltered VCF files `[ val(meta), path(*.vcf.gz) ]` *(if applicable)*
- `versions`: Channel with software versions `[ path(versions.yml) ]`

## Modules Used

1. **SPLIT_BED**: Splits BED file into equal chunks
2. **CALLING_VARDICT_CHUNK**: Runs VarDict on each chunk (parallel jobs)
3. **MERGE_VARDICT_RESULTS**: Merges all VCF and raw TSV chunk results into final VCF and a single gzipped raw TSV file

## Configuration

Configure via `ext.args` and `ext.filter_args` in `conf/modules.config`:

```nextflow
withName: '.*:BAM_CALL_VARDICT_PARALLEL:CALLING_VARDICT_CHUNK' {
    ext.args         = "-f 0.0 -r 1 -m 9999 -P 0 -p -z 1"
    ext.filter_args  = "-A -E -f 0.0 -p 10 -m 20 -v 1"
}
```

## Performance

For a 10GB BAM file with 48 chunks:
- **Sequential**: ~8 hours on 1 node with 48 CPUs
- **Parallel**: ~2 hours on 48 nodes with 4 CPUs each

## Example

```bash
nextflow run main.nf \
    --vardict_parallel_mode \
    --vardict_chunks 50 \
    --max_cpus 4
```


## Notes

- Only merged/final outputs (VCF and raw TSV.gz) are published to the output directory; intermediate chunk files are not published.
- The merged raw TSV file is created by concatenating all chunk raw TSVs and compressing to `.tsv.gz`.

This creates 50 independent jobs, each using up to 4 CPUs.
