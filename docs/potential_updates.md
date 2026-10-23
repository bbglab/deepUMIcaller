# Potential Updates

This file records options explored for improving individual pipeline steps that were
**not** implemented. Each entry documents the motivation, the alternatives considered,
and why the current implementation was kept (or what a future change would require).

## Replace `bamUtil trimBam` (TRIMBAM) with a faster, parallel tool

**Status: explored, not implemented — TRIMBAM retained.**

### Current implementation

- Module: `modules/nf-core/bamutil/trimbam/main.nf` (process `BAMUTIL_TRIMBAM`),
  container `biocontainers/bamutil:1.0.15--h2e03b76_1`.
- Applied in `workflows/deepumicaller.nf` only when `params.left_clip > 0` or
  `params.right_clip > 0`, on the **unmapped** BAM produced by UMI correction,
  before alignment.
- The module passes `-L $trim_left -R $trim_right` (from `ext.trim_left` /
  `ext.trim_right`, set in `conf/modules.config`) and **no `--clip` flag**, so
  trimBam runs in its default mode: trimmed bases are **masked to `N` with quality
  `!`** — read length is preserved, bases are not dropped.
- trimBam is single-threaded.

### Why "same output" is subtle

Downstream logic depends on the N-masking semantics:

- `modules/local/fgumi/filterconsensusreads/main.nf` computes
  `Ns_per_read = (params.maxN_per_read + params.left_clip + params.right_clip) /
  params.read_length_minus_tag` and passes it as `--max-no-call-fraction`. This
  budget assumes the trimmed bases are still present in the read as `N`s.
- Any replacement that **drops** bases (hard-clip, FASTQ-level trimming) changes
  read length and would require reworking this threshold, the
  `read_length_minus_tag` accounting, and validation that alignment and consensus
  calling are unaffected.

### Options explored

| Option | Output equivalence | Parallelism | Verdict |
| --- | --- | --- | --- |
| **A. Fold trimming into `fgumi extract`** via read structures (prepend an `S` segment, e.g. `10S+T` for a 10 bp left clip) | **No** — `S` ("Skip") drops the bases from `SEQ` entirely instead of N-masking; read length shrinks | Yes — `extract` is fully multi-threaded and already runs in the pipeline | Rejected as-is: conflicts with the `Ns_per_read` logic above. Would require updating that threshold and re-validating alignment/consensus. Attractive only if the N-masking semantics are deliberately abandoned. |
| **B. `fgumi clip --clipping-mode soft-with-mask`** with `--read-one-five-prime ${left_clip} --read-one-three-prime ${right_clip} --read-two-five-prime ${left_clip} --read-two-three-prime ${right_clip}` | **Yes (closest)** — soft-with-mask soft-clips *and* masks the clipped bases to `N` with minimum quality, matching trimBam's default mode | Yes — multi-threaded; the pipeline already uses `fgumi clip` post-consensus with the same left/right parameters (`conf/modules.config`, `CLIPBAM.*` `ext.extra_clipping`) | Best drop-in candidate. Caveats: requires `--reference` (a no-op for unmapped reads, but the argument must be supplied) and queryname-sorted/grouped input (the unmapped BAM from `fgumi extract` is query-grouped by construction). Would need a new local module plus a byte-level output comparison against trimBam on real data before switching. |
| **C. Custom parallel masking script** (e.g. `samtools view` + Python multiprocessing that sets trimmed bases to `N`/`!`) | **Yes** — semantics can be made exactly identical | Yes — bounded by chosen worker count | Fallback if Option B's output diverges. Cost: new code to write and maintain, no upstream bug fixes, and it duplicates what Option B already provides. |
| **D. FASTQ-level trimming before `fgumi extract`** (e.g. fastp `--trim_front1`, seqkit) | **No** — these drop bases rather than mask them | Yes | Rejected: same read-length/`Ns_per_read` conflict as Option A, and it adds a process rather than removing one. |
| **E. Other BAM-level clippers** (fgbio `ClipBam`, jvarkit, bamtools) | fgbio `ClipBam` can mask but is Java/single-threaded like trimBam; jvarkit/bamtools have no maintained N-masking trim equivalent | Mostly no | Rejected: no speed or parallelism win over the status quo. |

### Recommendation if revisited

Option B (`fgumi clip --clipping-mode soft-with-mask`) is the only candidate that
satisfies all three constraints (same output, faster, parallel) without touching
downstream thresholds. Suggested validation plan before switching:

1. Add a local `FGUMI_CLIPBAM_UNMAPPED`-style module fed by `ch_ref_fasta` and the
   corrected unmapped BAM, with the same `left_clip`/`right_clip` parameters.
2. Run both TRIMBAM and the new module on a test sample and diff the BAMs
   record-by-record (SEQ, QUAL, flags, read length) — expect identity apart from
   `@PG` header lines.
3. Confirm the `filterconsensusreads` `Ns_per_read` threshold still holds, then
   remove the `bamutil` dependency from `modules.json` and the nf-core module tree.

Note: fgumi is already pinned at 0.7.0 (containers) / 0.8.0 (conda) in this
pipeline, so Option B adds no new dependency.
