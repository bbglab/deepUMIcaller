# bbglab/deepUMIcaller: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Migrate UMI tooling from `fgbio` to `fgumi` (Fulcrum Genomics UMI toolkit): `fgumi fastq`, `fgumi zipper`, `fgumi sort`, `fgumi clip` and `fgumi group` replace the corresponding `fgbio` subcommands, with container updates through `fgumi 0.7.0` ([#209]).
- Switch the aligner in `ALIGN_BAM` from `bwa` to `bwa-mem3` (`bwa-mem3 mem`).
- New `preprocessed_crams` entrypoint that accepts pre-aligned CRAMs and retags them with `fgumi retag` ([#209]).
- Compute UMI collisions (`compute_umi_collisions.py`) and report them in the QC outputs.
- Trim reads containing `N` bases before alignment ([#202]).
- Support multiple pre-built reference indices ([#204]).

### Changed

- Tune memory allocation for the heavy BAM-processing steps.
- General linting fixes across modules and configs.

### Fixed

- Fix a bug in the UMI collisions computation.

## v1.1.3 - 2026-07-31

### Added

- Internal FASTQ splitting with `seqkit split2` ([#196]).

### Fixed

- Depth correction applied to the mutations-per-position ratios.

## v1.1.2 - 2026-06-03

### Added

- Expand internal panel ([#188]).
- Report non-mutated reads ([#190]).
- Conditional `pileupbams` output ([#187]).

### Changed

- Documentation update ([#193]).

## v1.1.1 - 2026-04-09

### Changed

- Rework the output directory structure ([#179], [#183]).
- Schema now requires integer sample IDs.
- Auto-detect split inputs.

### Fixed

- Handle silent VarDict failures.
- Update test data and docs ([#135]).

## v1.1.0 - 2026-02-06

### Added

- DeepClone support.
- VarDict true parallelization with `vardict_chunks` ([#161]).
- Unmapped BAM step ([#159]).
- BED-related masks ([#155]).
- ClipBAM for the all-molecules BAM ([#176]).
- Mutated reads output ([#173]).
- Sample validation ([#175]).
- Early MultiQC duplex report.
- Low-complexity BED mask generation script.

### Changed

- `as_minus_xs` rewritten as a Click CLI with a two-pass optimization.
- Base config refactor ([#174]).
- Default read size set to 150.

### Fixed

- `parent_dna` merging fixes ([#165]).

## v1.0.0 - 2025-10-08

### Added

- `nanoseq` BED-based mutation filters ([#141]).
- GPLv3 license ([#140]).
- Single confidence level family metrics ([#147]).
- deepCSA input generation.
- `parent_dna` run mode.

### Fixed

- Calling parallelization bug ([#153]).

## v0.5.0 - 2025-08-26

### Changed

- Documentation updates ([#128]).

## v0.4.0 - 2025-07-22

### Added

- Purine/pyrimidine split plotting ([#106]).
- Discarded BAM coverage ([#101]).

### Changed

- Reposition AM filtering ([#97]).
- BED file handling ([#95]).

## v0.3.0 - 2025-03-27

### Added

- AM filtering ([#90]).
- Custom read size support.

### Changed

- Refactor `n_rich` filtering ([#93]).
- VarDict parameter tuning.
- Optimize depth recomputation ([#87]).

## v0.2.0 - 2025-01-13

### Added

- Light mode ([#86]).

### Changed

- Samtools 1.20.
- 24.04 resource format.

### Fixed

- QC turn-off fixes.
- BED-to-BAM pairing fix.
- VCF reading bug fixes.

## v0.1.0-beta - 2024-01-25

### Added

- First fork release from nf-core/fastquorum, including the `callingbed` step ([#48]).
