---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# Codebase Concerns

**Analysis Date:** 2026-10-10

## Tech Debt

**Temporary fixes config promoted to permanent:**
- Issue: `conf/tmp_quick_fixes.config` sets `errorStrategy = 'ignore'` with `maxRetries = 1` for all processes matching `.*COVERAGE.*`. Failed coverage steps are silently swallowed and the pipeline continues, producing incomplete downstream results without any failure signal.
- Files: `conf/tmp_quick_fixes.config`
- Impact: Silent data loss — a failed COVERAGE process yields missing/partial outputs that later stages may treat as valid empties.
- Fix approach: Identify the root cause of the COVERAGE failures, apply a targeted fix, then either delete this config or scope the ignore to a specific process with a logged warning.

**Fork residue from nf-core/fastquorum:**
- Issue: Pipeline is a fork of nf-core/fastquorum but branding/assets were not fully migrated. `subworkflows/local/utils_nfcore_deepumicaller/main.nf` has unfilled citation/Zenodo placeholders (`TODO nf-core: Add Zenodo DOI for pipeline after first release`, lines 211-213).
- Files: `subworkflows/local/utils_nfcore_deepumicaller/main.nf`
- Impact: Missing citation/DOI metadata.
- Fix approach: Fill citation text, register a Zenodo DOI.
- **RESOLVED (2026-10-21):** The sendmail-template fork residue is fixed — `assets/sendmail_template.txt` now embeds `assets/bbglab_logo.png` (copied from `docs/images/`), and the previously missing `assets/email_template.txt`/`.html` were created (deepUMIcaller-branded). Email notifications are now fully wired: `main.nf` invokes `PIPELINE_COMPLETION` with the MultiQC report. The Zenodo DOI is also filled: `nextflow.config` manifest now declares `doi = '10.5281/zenodo.19487799'` (which also removes the `nodoi_text` hint from MultiQC methods text), the README badge is active, and `citation()` in `subworkflows/local/utils_nfcore_deepumicaller/main.nf` lists it.

**Stale CHANGELOG:**
- Issue: `CHANGELOG.md` still shows `## v1.0dev - [date]` with a literal `[date]` placeholder, while `nextflow.config` manifest declares `version = 'v1.1.0 Montseny'` (line 288). Release history is not tracked.
- Files: `CHANGELOG.md`, `nextflow.config`
- Impact: No traceability of what changed between releases; version mismatch between docs and manifest.
- Fix approach: Backfill changelog entries for v1.1.0 and adopt a release checklist step.

**Unimplemented MNV depth recomputation:**
- Issue: `bin/recompute_depth.py` line 276 passes an empty list with comment `// FIXME implement mutated reads for MNVs`; MNV records get `alt_dp_count = -1` placeholder values.
- Files: `bin/recompute_depth.py`
- Impact: MNV variants have incorrect/placeholder depth annotations in final outputs.
- Fix approach: Implement per-base read matching for multi-nucleotide variants, or explicitly drop/flag MNV records.

**Undocumented local modules:**
- Issue: Several local modules carry `TODO` markers for missing documentation and proper nf-core structure: `modules/local/filter/nrich/main.nf` (line 25), `modules/local/filter/nrich/meta.yml` (line 1), `modules/local/filtervcf/main.nf` (line 4), `modules/local/filter_reads/samtools/main.nf` (line 33 — ad-hoc samtools view instead of a formal module), `modules/local/sigprofiler/matrixgenerator/meta.yml` (template TODOs throughout), `subworkflows/local/postprocess_mutations/meta.yml`.
- Files: `modules/local/filter/nrich/main.nf`, `modules/local/filter/nrich/meta.yml`, `modules/local/filtervcf/main.nf`, `modules/local/filter_reads/samtools/main.nf`, `modules/local/sigprofiler/matrixgenerator/meta.yml`, `subworkflows/local/postprocess_mutations/meta.yml`
- Impact: Modules are hard to reason about; meta.yml stubs break nf-core lint expectations.
- Fix approach: Fill in meta.yml descriptions and refactor the ad-hoc samtools filter into a named, configurable module.

**fgumi migration incomplete (in-flight branch):**
- Issue: Current branch `replace-fgbio-with-fgumi-bwamem3` migrates fgbio → fgumi, but remnants remain: `modules/local/fgumi/clipbam/main.nf` line 23 has `// FIXME : this should be switched to fgumi sort` (currently pipes `fgumi sort | fgumi clip` through /dev/stdin), and `tests/test_data/expected_output/fgbio/` still holds legacy expected outputs alongside new top-level VCFs.
- Files: `modules/local/fgumi/clipbam/main.nf`, `tests/test_data/expected_output/fgbio/`, `modules.json`
- Impact: Two-generation tooling coexists; the stdin pipe in CLIPBAM adds an unnecessary failure point and complicates resource attribution.
- Fix approach: Complete the fgumi sort switch, remove or archive the `fgbio/` expected-output directory, and merge the branch before further divergence.

## Known Bugs

**`bool(sys.argv[6])` misparses split flag in filtervcf:** — RESOLVED 2026-10-21
- Symptoms: Any non-empty string argument — including `"false"` or `"0"` — evaluates to `True`, so `split_type` is enabled unintentionally when a 6th argument is passed.
- Files: `bin/filtervcf.py` (lines 227-229)
- Trigger: Invoking `filtervcf.py` with 6 positional args where the last is any string (e.g., `"false"`).
- Workaround: Omit the 6th argument entirely (the `except` branch defaults to `False`).
- Resolution: Migrated to argparse (`--split-type` store_true flag); splitting is now only enabled by the explicit flag. Call sites in `modules/local/filtervcf/main.nf` updated.

**Bare `except:` clauses swallow all errors:** — PARTIALLY RESOLVED 2026-10-21
- Symptoms: `bin/family_size_plots.py` lines 91 and 154 use bare `except:` — a failed `np.argmax` silently sets `peak_size = 0`, and unreadable duplex metrics files are skipped with only a print (no non-zero exit). `bin/filtervcf.py` line 228 uses bare `except:` for argument parsing.
- Files: `bin/family_size_plots.py`, `bin/filtervcf.py`
- Trigger: Empty/NaN plot data; corrupt or missing metrics TSVs.
- Workaround: None — outputs are silently degraded.
- Resolution: The `filtervcf.py` bare `except:` (argument parsing) is gone after the argparse migration. The `family_size_plots.py` bare `except:` clauses remain open.

**Broad `except Exception` blocks hide data problems:**
- Symptoms: `bin/summarize_muts_per_cycle.py` has five broad exception handlers (lines 63, 85, 103, 120, 165) that skip malformed rows/files without failing.
- Files: `bin/summarize_muts_per_cycle.py`
- Trigger: Malformed mpileup or metrics input.
- Workaround: None — summaries silently omit data.

## Security Considerations

**Institutional infrastructure paths committed to the repo:**
- Risk: `conf/general_files_IRB.config` and `conf/test.config` hardcode internal cluster paths (`/data/bbg/datasets/...`, singularity `cacheDir`/`libraryDir`), exposing internal infrastructure layout in a public-facing repo.
- Files: `conf/general_files_IRB.config`, `conf/test.config`
- Current mitigation: Paths are data locations, not credentials; no secrets were found in any config (checked `conf/*.config`, `tower.yml`).
- Recommendations: Move site-specific paths to a private/institutional config profile (nf-core `custom_config` pattern, already partially wired via `params.custom_config_base` in `nextflow.config` lines 124-125) and keep only placeholders in-repo.

**No CI pipeline:**
- Risk: `.github/workflows/` does not exist — no automated linting (nf-core lint, prettier), no nf-test execution on push/PR. All validation is manual on the IRB cluster.
- Files: (missing `.github/workflows/`)
- Current mitigation: Manual nf-test runs (`nf-test.config`, `tests/main.nf.test`).
- Recommendations: Add nf-core GitHub Actions workflows (lint, test) adapted for the fork.

**Email notification template with no parameter wiring:**
- Risk: `assets/sendmail_template.txt` exists but no `sendmail`/notification parameters were found in `nextflow.config` or workflow files — dead template referencing fork-origin logos.
- Files: `assets/sendmail_template.txt`
- Current mitigation: Unused code path.
- Recommendations: Either wire up notification params or remove the template.
- **RESOLVED (2026-10-21):** Notifications are now wired, not removed. `main.nf` invokes `PIPELINE_COMPLETION` (email + summary) after `DEEPUMICALLER`, passing `params.email`/`params.email_on_fail`/`params.plaintext_email` and the MultiQC report (new `emit:` on `DEEPUMICALLER`). Created the missing `assets/email_template.txt`/`.html` and re-pointed `assets/sendmail_template.txt` at `assets/bbglab_logo.png`. `hook_url` (webhook) intentionally left unwired (passed as `null`). Delivery is best-effort and fail-safe: sendmail exit status is checked, hung MTAs are killed after 60 s, and any delivery error only logs a warning — the run result is never affected. Caveat: actual delivery still depends on the host MTA relaying external mail — verify with a real `--email` run.

## Performance Bottlenecks

**Silenced COVERAGE failures mask resource under-provisioning:**
- Problem: The `errorStrategy = 'ignore'` patch in `conf/tmp_quick_fixes.config` exists because COVERAGE processes fail — likely memory/time limits rather than data issues.
- Files: `conf/tmp_quick_fixes.config`, `conf/base.config` (label definitions)
- Cause: Under-sized label resources for coverage computation on large panels.
- Improvement path: Profile the failing COVERAGE tasks, bump the relevant label memory, remove the ignore.

**fgumi sort/clip pipeline through /dev/stdin:**
- Problem: `modules/local/fgumi/clipbam/main.nf` chains `fgumi sort --output - | fgumi clip --input /dev/stdin`, serializing two heavy steps in one task and doubling peak memory within a single allocation.
- Files: `modules/local/fgumi/clipbam/main.nf`
- Cause: FIXME acknowledges this should be a single `fgumi sort` invocation.
- Improvement path: Switch to the integrated fgumi sort (per the in-file FIXME) or split into two processes with independent retry.

**Chromosome splitting re-reads the full BAM per chromosome:**
- Problem: `modules/local/splitbamchrom/main.nf` runs `samtools view -c` on the whole BAM plus per-chromosome views and stats — roughly 4+ full-BAM passes per sample.
- Files: `modules/local/splitbamchrom/main.nf`
- Cause: Stats collection interleaved with splitting in one shell loop.
- Improvement path: Compute idxstats once, derive counts from it, and split in a single pass.

## Fragile Areas

**Main workflow branching (`workflows/deepumicaller.nf`, 652 lines):**
- Files: `workflows/deepumicaller.nf`
- Why fragile: Dense conditional logic on `params.split_by_chrom` (lines 310, 361, 379, 421) and `meta.parent_dna` grouping (lines 463-490) with inline meta-map mutation; channel topology changes for one mode can silently break the other.
- Safe modification: Change one branch at a time; run both `split_by_chrom=true` and `false` nf-test cases (`tests/main.nf.test` TEST 1 and TEST 2 families) after any edit.
- Test coverage: Pipeline-level nf-tests exist for both modes; no unit tests on the meta-manipulation logic itself.

**Chunk-ID regex in VarDict chunking:**
- Files: `modules/local/calling_vardict/calling_vardict_chunk/main.nf`
- Why fragile: `meta.id.replaceAll(/_chunk_[^_]+$/, '')` reconstructs the original sample ID by string surgery; a sample ID legitimately containing `_chunk_` would be corrupted. Inline `awk` fallback when `params.use_teststrandbias` is false duplicates VarDict's filtering logic in shell.
- Safe modification: Pass the original sample ID explicitly through meta instead of regex-stripping; keep the awk path in sync with `teststrandbias.R` behavior when touching either.
- Test coverage: Covered only indirectly via pipeline tests.

**Species-hardcoded chromosome lists:**
- Files: `modules/local/splitbamchrom/main.nf`
- Why fragile: Chromosome lists are hardcoded per `params.vep_species` (`homo_sapiens`, `mus_musculus`) with a hard `error` for anything else; contig naming assumptions (`^chr${chrom}$|^chr${chrom}_`) break on non-chr-prefixed or non-standard assemblies.
- Safe modification: Derive chromosomes from the BAM idxstats or a params-supplied list rather than species lookup.
- Test coverage: Only human test data exercised.

**Positional-argument Python entry points:** — RESOLVED 2026-10-21
- Files: `bin/filtervcf.py` (6 positional args), `bin/VCF_Parser.py`, `bin/add_filter_from_bed.py`
- Why fragile: Scripts parse `sys.argv` positionally with no argparse validation in some cases; argument order changes break callers silently (see the `bool(sys.argv[6])` bug above).
- Safe modification: Migrate to argparse with named flags (pattern already used by `bin/check_samplesheet.py` `parse_args`).
- Test coverage: None (no Python unit tests).
- Resolution: `bin/filtervcf.py` and `bin/add_filter_from_bed.py` migrated to argparse with named flags (`parse_args` pattern); Nextflow call sites updated (`modules/local/filtervcf/main.nf`, `modules/local/filter/from_bed/main.nf`). `bin/VCF_Parser.py` needs no migration — it is a pure library (no `__main__`, no `sys.argv`), imported by `bin/count_muts_per_cycle.py`. Note `bin/compare_vcfs.py` still uses raw `sys.argv` (see CONVENTIONS.md Style 3) but was out of scope here.

## Scaling Limits

**Resource ceilings tuned to one cluster:**
- Current capacity: `nextflow.config` sets `max_memory = '1024.GB'`, `max_cpus = 192` (lines 133-134); `conf/local.config` caps at 16 cpus / 60.GB / 240.h.
- Limit: Portability — running elsewhere requires overriding these plus the `/data/bbg` paths.
- Scaling path: Use nf-core institutional config profiles instead of in-repo hardcoded ceilings.

**`qc_memory_intensive` label at 170.GB:**
- Current capacity: `conf/base.config` line 47-49 allocates 170.GB for QC-memory-intensive tasks.
- Limit: Nodes with less RAM will queue indefinitely or fail.
- Scaling path: Verify which processes use this label and whether the allocation reflects measured peaks.

## Dependencies at Risk

**nf-core modules pinned to floating `master` branch SHAs:**
- Risk: `modules.json` pins every nf-core module to a `git_sha` on `branch: master` with no nf-core release tag; `nf-core modules update` can pull breaking changes, and SHAs drift from tested combinations.
- Files: `modules.json`, `modules/nf-core/*`, `subworkflows/nf-core/*`
- Impact: Reproducibility depends on GitHub retaining those exact SHAs; updates are all-or-nothing per module.
- Migration plan: Adopt nf-core module versioning discipline (`nf-core modules update` with review, or pin to a released nf-core modules tag).

**Mixed container sources for local modules:**
- Risk: `modules/local/filtervcf/main.nf` pins `biocontainers/seaborn:0.12.2_cv1` (a plotting library container) while relying on `bin/filtervcf.py` mounted from the pipeline `bin/` — the container does not contain the tool it runs. Conda fallback `anaconda::seaborn=0.12.2` differs from the container environment.
- Files: `modules/local/filtervcf/main.nf`, `bin/filtervcf.py`
- Impact: Container/conda divergence; the in-file TODO notes the proper container is not yet in biocontainers.
- Migration plan: Build a dedicated container for filtervcf dependencies when available.

**fgbio/fgumi transition window:**
- Risk: Until branch `replace-fgbio-with-fgumi-bwamem3` merges, consensus-calling tooling is in flux; `modules.json` still pins `fgbio/*` modules while local `fgumi/*` modules replace them.
- Files: `modules.json`, `modules/nf-core/fgbio/` (if present), `modules/local/fgumi/*`
- Impact: Divergent consensus outputs between branches; expected outputs must be regenerated per tool generation.
- Migration plan: Merge the migration branch, regenerate `tests/test_data/expected_output/`, remove legacy `fgbio/` fixtures.

## Missing Critical Features

**No CI automation:**
- Problem: No `.github/workflows/` — linting and tests run only manually.
- Blocks: Confident merges; the uncommitted `conf/test.config` and modified `tests/test_data/expected_output/Bsample.vcf` in the working tree show manual-test-driven drift.

**MNV support in depth recomputation:**
- Problem: `bin/recompute_depth.py` FIXME (line 276) — MNV reads are not counted.
- Blocks: Accurate depth annotation for MNV calls in final VCFs.

**Pipeline DOI/citation:**
- Problem: Zenodo DOI placeholder unfilled in `subworkflows/local/utils_nfcore_deepumicaller/main.nf` (lines 211-213).
- Blocks: Proper citation of released versions.
- **RESOLVED (2026-10-21):** DOI `10.5281/zenodo.19487799` set in the `nextflow.config` manifest, README badge, `citation()` helper, and `CITATIONS.md`.

## Test Coverage Gaps

**Zero unit tests for `bin/` Python scripts:**
- What's not tested: All 3,187 lines across 12 scripts (`bin/recompute_depth.py` 682 lines, `bin/count_muts_per_cycle.py` 448, `bin/check_samplesheet.py` 318, `bin/compute_umi_collisions.py` 289, `bin/filtervcf.py` 271, `bin/VCF_Parser.py` 258, etc.). No `test_*.py` files exist anywhere in the repo. `bin/test_utilities/` holds ad-hoc comparison helpers (`compare_vcfs_detailed.py`, `validate_vcf_biology.py`), not a test suite.
- Files: `bin/*.py`
- Risk: The `bool(sys.argv[6])` bug and bare-except swallowing would be caught trivially by unit tests; logic regressions surface only as full-pipeline VCF diffs.
- Priority: High — start with `filtervcf.py`, `recompute_depth.py`, `check_samplesheet.py`.

**nf-test suite is single-file and pipeline-level only:**
- What's not tested: Individual local modules/subworkflows in isolation; `nf-test.config` explicitly ignores `modules/nf-core/**` and `subworkflows/nf-core/**` (acceptable — upstream CI), but there are no per-module tests for the ~20 `modules/local/*` processes either.
- Files: `tests/main.nf.test` (863 lines, 11 tests), `nf-test.config`
- Risk: A change to one local module requires a full pipeline run to validate.
- Priority: Medium — add module-level nf-tests for the fragile areas listed above.

**Test result artifacts uncommitted and expected outputs in flux:**
- What's not tested: `tests/2026-10-07.results.csv` and `tests/2026-10-09.results.csv` are untracked; `tests/test_data/expected_output/Bsample.vcf` is modified and `tests/test_data/expected_output/fgbio/Bsample.vcf` is untracked in the working tree.
- Files: `tests/`, `tests/test_data/expected_output/`
- Risk: Expected outputs and reality can drift without a commit gate; unclear which expected-output generation is canonical during the fgbio→fgumi migration.
- Priority: Medium — commit or clean these before the migration branch merges.

---

*Concerns audit: 2026-10-10*
