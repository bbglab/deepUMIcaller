# bbglab/deepUMIcaller: Output

## Introduction

This document collects specific details and comments on tools used in the pipeline or decisions made along the process.

<!-- ## Document overview -->

<!-- - [Basic mutations QC](#basic-mutations-qc)
- [Interal processing files](#internal-files) -->

## VarDict

VarDict is not calling any variants in the first and last 5 bps of a read.

This should be taken into account when setting the clipping parameters, but also note that since we do the pileup for redefining the frequency of each mutation (ALT_DEPTH) the mutations present in these first positions would still be taken into account.

## Read-end clipping (TrimBam)

When `left_clip` and/or `right_clip` are set, deepUMIcaller trims read ends using [`bamutil trimBam`](https://genome.sph.umich.edu/wiki/BamUtil:_trimBam) with the `--clip` option. This clips the ends (soft-clipping) instead of replacing trimmed bases with `N`, so the clipped positions are removed from downstream read content and do not inflate the `N` fraction used by `fgbio FilterConsensusReads`.

## Internal panel expansion (EXPAND_PANEL)

The first process executed in every run pads the user-supplied `--targetsfile` BED by ±250 bp on each side. The expanded BED is then used for downstream coverage/capture/calling steps, including the on-target version of `fgbio CollectDuplexSeqMetrics`. Users should provide their raw targeted regions; pre-extending the BED is not necessary and would result in double-padding.

## Two-pass depth recomputation (PATCHDP / PATCHDPALL)

VarDict's depth fields are unreliable for duplex consensus data. After calling, deepUMIcaller replaces them with counts obtained from `samtools mpileup`. A first pass is run over the duplex BAM (`PATCHDP`); a second pass is run over the all-molecules BAM (`PATCHDPALL`) to add depth/VAF information for non-duplex reads. The final VCF therefore carries information to compute four VAF values per mutation: `VAF_vd` (VarDict), `VAF` (duplex), `VAF_AM` (all-molecules) and `VAF_ND` (no-duplex).

## fgbio FilterConsensusReads — two stringencies

The pipeline runs `fgbio FilterConsensusReads` twice on the same all-molecules BAM:

- **All-molecules filter** (`FILTERCONSENSUSREADSAM`, default `--min-reads 2 1 0`, no strand agreement required) — keeps both duplex and single-strand consensus reads (with at least 2 copies in a single strand); produces the BAM used to compute `VAF_AM` and `VAF_ND` at postprocessing time.
- **Duplex filter** (`FILTERCONSENSUSREADSDUPLEX`, default `--min-reads 4 2 2`, `--require-single-strand-agreement true`) — keeps only high-quality duplex consensus reads; produces the BAM used for variant calling.
