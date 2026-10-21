---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# Codebase Structure

**Analysis Date:** 2026-10-10

## Directory Layout

```
deepUMIcaller/
├── main.nf                    # Entry point: PIPELINE_INITIALISATION + DEEPUMICALLER
├── nextflow.config            # Global params, profiles, plugins, manifest
├── nextflow_schema.json       # Param schema for nf-schema validation
├── modules.json               # nf-core modules dependency manifest (nf-core/tools)
├── nf-test.config             # nf-test runner config
├── tower.yml                  # Seqera Tower launch config
├── workflows/
│   └── deepumicaller.nf       # Main workflow: all channel wiring (652 lines)
├── subworkflows/
│   ├── local/                 # Pipeline-specific subworkflows
│   │   ├── input_check.nf                 # Samplesheet → meta channels
│   │   ├── bam_call_vardict_parallel/     # Chunked parallel VarDict
│   │   ├── postprocess_mutations/         # Pileup/patch/filter/muts-per-pos
│   │   ├── vcf_annotate_all/              # VEP annotation wrapper
│   │   └── utils_nfcore_deepumicaller/    # Init + completion + logo/citation
│   └── nf-core/               # Vendored nf-core subworkflows
│       ├── utils_nextflow_pipeline/
│       ├── utils_nfcore_pipeline/
│       ├── utils_nfvalidation_plugin/
│       └── vcf_annotate_ensemblvep/
├── modules/
│   ├── local/                 # Custom processes (see Directory Purposes)
│   │   ├── utils.nf           # Shared Groovy functions (process_bams, clean_chr_names)
│   │   └── <tool>/<process>/main.nf (+ meta.yml)
│   └── nf-core/               # Vendored nf-core modules
│       ├── bamutil/ bedtools/ custom/ ensemblvep/ fastp/ fastqc/
│       ├── fgumi/ (groupreadsbyumi, callduplexconsensusreads)
│       ├── multiqc/ picard/ qualimap/ samtools/
├── conf/                      # Config profiles (see Directory Purposes)
├── bin/                       # Python helpers, auto-added to PATH
├── assets/                    # Samplesheet schema, MultiQC config, templates
├── docs/                      # usage.md, output.md, input_scenarios.md, tools.md
├── tests/                     # nf-test pipeline tests + fixtures
│   ├── main.nf.test
│   ├── nextflow.config
│   └── test_data/{input,expected_output}/
├── testing/                   # (empty)
├── work/                      # Nextflow work dir (generated, uncommitted)
└── .planning/                 # GSD planning docs (this directory)
```

## Directory Purposes

**`workflows/`:**
- Purpose: top-level workflow definitions
- Contains: `deepumicaller.nf` — the single orchestration workflow with all imports and `params.step` gates
- Key files: `workflows/deepumicaller.nf`

**`subworkflows/local/`:**
- Purpose: pipeline-specific multi-process groups
- Contains: `input_check.nf` (flat file), `bam_call_vardict_parallel/main.nf`, `postprocess_mutations/main.nf`, `vcf_annotate_all/main.nf`, `utils_nfcore_deepumicaller/main.nf`
- Key files: `subworkflows/local/postprocess_mutations/main.nf` (largest, 227 lines)

**`subworkflows/nf-core/`:**
- Purpose: vendored nf-core subworkflows (copied, not fetched via modules.json at runtime)
- Contains: pipeline utils (init/completion/validation) and VEP annotation
- Key files: `subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf`

**`modules/local/`:**
- Purpose: custom process definitions, one directory per tool
- Contains:
  - `align_bam/` — `ALIGN_BAM` (fgumi fastq | bwa-mem3 | fgumi zipper | fgumi sort)
  - `fgumi/` — `fastqtobam/`, `correctumis/`, `retag/`, `clipbam/`, `filterconsensusreads/`, `collectduplexseqmetrics/`, `sort/`
  - `calling_vardict/` — `split_bed/`, `calling_vardict_chunk/`, `merge_vardict_results/`
  - `filter/` — `from_bed/` (FILTER_FROM_BED), `nrich/` (FILTER_N_RICH)
  - `filter_reads/` — `asminusxs/`, `samtools/`
  - `filtervcf/`, `count_ns/`, `patchdepth/`, `tabix_mpileup/`, `mutsperpos/{compute,summarize}/`
  - `familymetrics/`, `umicollisions/`, `createbed/`, `expand_regions/`
  - `mergebam/`, `splitbamchrom/`, `splitfastq/`, `unmap/`
  - `sigprofiler/matrixgenerator/`, `bedtools/merge/`
  - `samplesheet_check.nf` (flat file, not a directory)
  - `utils.nf` — shared Groovy functions
- Key files: `modules/local/align_bam/main.nf`, `modules/local/utils.nf`

**`modules/nf-core/`:**
- Purpose: vendored nf-core modules (fastp, fastqc, samtools, picard, qualimap, bedtools, multiqc, ensemblvep, bamutil, custom/dumpsoftwareversions, fgumi groupreadsbyumi + callduplexconsensusreads)
- Contains: standard nf-core `main.nf` + `meta.yml` per process
- Key files: `modules/nf-core/fgumi/groupreadsbyumi/main.nf`, `modules/nf-core/samtools/mpileup/main.nf`

**`conf/`:**
- Purpose: layered Nextflow configuration
- Contains:
  - `base.config` — default resources, labels, retry strategy (loaded always)
  - `modules.config` — per-process `ext.args`, default publishDir, disabled publishing for intermediates
  - `results_outputs.config` — publishDir paths per process (metrics/, mutations_vcf/, processing_files/)
  - `exome.config` — exome-tuned resource labels
  - `local.config` — local executor limits (16 cpus / 60 GB)
  - `igenomes.config` — iGenomes reference paths (GRCh37/GRCh38)
  - `test.config`, `test_full.config` — test profiles
  - `general_files_IRB.config` — IRB Barcelona cluster paths + singularity cache (profile `irbcluster`)
  - `tmp_quick_fixes.config` — temporary process patches
- Key files: `conf/base.config`, `conf/modules.config`, `conf/results_outputs.config`

**`bin/`:**
- Purpose: Python helper scripts; Nextflow puts `bin/` on PATH inside every process
- Contains: `filtervcf.py`, `VCF_Parser.py`, `compute_umi_collisions.py`, `count_muts_per_cycle.py`, `summarize_muts_per_cycle.py`, `family_size_plots.py`, `as_minus_xs.py`, `add_filter_from_bed.py`, `add_filter_nrich.py`, `recompute_depth.py`, `check_samplesheet.py`, `compare_vcfs.py`, `test_utilities/`
- Key files: `bin/filtervcf.py`, `bin/compute_umi_collisions.py`

**`assets/`:**
- Purpose: static pipeline assets
- Contains: `schema_input.json` (samplesheet schema), `multiqc_config.yml`, `samplesheet.csv` (example), `email_template.txt`/`email_template.html` (completion email body), `sendmail_template.txt` (MIME wrapper, embeds `bbglab_logo.png`), `bbglab_logo.png`, `generate_low_complex_rep_bed.py`, `images/`
- Key files: `assets/schema_input.json`

**`tests/`:**
- Purpose: nf-test pipeline-level functional tests
- Contains: `main.nf.test` (multi-scenario suite), `nextflow.config`, `test_data/input/*.csv` (per-step samplesheets), `test_data/expected_output/*.vcf`, dated results CSVs (`2026-10-07.results.csv`, `2026-10-09.results.csv`)
- Key files: `tests/main.nf.test`

**`docs/`:**
- Purpose: user documentation
- Contains: `usage.md`, `output.md`, `input_scenarios.md`, `test_data.md`, `tools.md`, `README.md`, `images/`

## Key File Locations

**Entry Points:**
- `main.nf`: pipeline entry; declared as `mainScript` in `nextflow.config:247`
- `workflows/deepumicaller.nf`: `DEEPUMICALLER` workflow definition
- `tests/main.nf.test`: nf-test suite entry

**Configuration:**
- `nextflow.config`: global params, profiles (docker/singularity/conda/test/local/exome/irbcluster...), manifest, plugins
- `nextflow_schema.json`: parameter schema (nf-schema validation + help)
- `conf/base.config`: resource labels and retry policy
- `conf/modules.config`: per-module `ext.args` and default publishing
- `conf/results_outputs.config`: output publishing paths
- `assets/schema_input.json`: samplesheet column schema

**Core Logic:**
- `workflows/deepumicaller.nf`: all channel topology and step gating
- `subworkflows/local/input_check.nf`: meta-map construction (`create_input_channel`)
- `subworkflows/local/bam_call_vardict_parallel/main.nf`: parallel calling logic
- `subworkflows/local/postprocess_mutations/main.nf`: filter chain logic
- `modules/local/utils.nf`: shared channel helper functions

**Testing:**
- `tests/main.nf.test`: pipeline tests (basic, multi-file, multi-sample, split_by_chrom, per-step entries)
- `tests/test_data/input/`: per-step samplesheet fixtures
- `tests/test_data/expected_output/`: reference VCFs for ≥99% precision comparison
- `bin/test_utilities/`: VCF comparison/validation helpers (`compare_vcfs_detailed.py`, `validate_vcf_biology.py`)

## Naming Conventions

**Files:**
- Module/subworkflow logic: `main.nf` inside `<tool>/<process>/` directories
- Metadata: `meta.yml` alongside each vendored module's `main.nf`
- Flat single-process modules: snake_case `.nf` (`modules/local/samplesheet_check.nf`, `subworkflows/local/input_check.nf`)
- Python helpers: snake_case (`bin/compute_umi_collisions.py`)
- Configs: lowercase with `.config` extension in `conf/`

**Directories:**
- `modules/<local|nf-core>/<tool>/<process>/` — tool grouped, process nested
- Multi-process tools nest one level deeper: `modules/local/fgumi/clipbam/`, `modules/local/calling_vardict/split_bed/`

**Processes:**
- UPPER_SNAKE_CASE, often tool-prefixed: `FGUMI_FASTQTOBAM`, `FGUMI_CLIPBAM`, `SAMTOOLS_MPILEUP`, `BEDTOOLS_MERGE`, `PICARD_MERGESAMFILES`
- Aliased at import to role-based names: `ALIGN_BAM as ALIGNRAWBAM`, `FGUMI_CLIPBAM as CLIPBAMAM`, `SAMTOOLS_MPILEUP as PILEUPBAMALL`
- Labels (resource classes): kebab-ish lowercase snake: `alignment_intensive`, `variant_calling`, `qc_memory_intensive`

**Channels/meta:**
- Channel vars: `ch_` prefix (`ch_multiqc_files`, `ch_targetsfile`) or descriptive snake_case (`aligned_raw_bam`, `cons_duplex_bam`)
- Meta keys: `id`, `sample`, `parent_dna`, `read_structure`, `umi_file`, `chunk_name`

## Where to Add New Code

**New Process (custom tool):**
- Implementation: `modules/local/<tool>/<process>/main.nf` with `process <TOOL>_<PROCESS> { ... }`, `tag "$meta.id"`, `label '<resource>'`, `conda`/`container` directives, `tuple val(meta), path(...)` I/O, `versions.yml` with `topic: versions`, optional `stub:`
- Wire it: import in `workflows/deepumicaller.nf` (alias if reused), call inside the appropriate `if (params.step in [...])` block
- Resources: add/extend a `withLabel:` block in `conf/base.config` (and `conf/exome.config` if exome-relevant)
- Publishing: add a `withName:` block in `conf/results_outputs.config` (or `conf/modules.config` for args/disabled publishing)

**New Subworkflow:**
- Implementation: `subworkflows/local/<name>/main.nf` with `workflow <NAME> { take: ... main: ... emit: ... }`
- Wire it: `include { <NAME> } from '../subworkflows/local/<name>/main'` in `workflows/deepumicaller.nf`

**New Python Helper:**
- Implementation: `bin/<snake_case>.py` (auto on PATH; add `click` if CLI-style, matching `bin/recompute_depth.py`, or the argparse `parse_args(argv=None)` pattern matching `bin/filtervcf.py`/`bin/add_filter_from_bed.py`)
- No registration needed; call by name in a process `script:` block

**New Parameter:**
- Default: `params` block in `nextflow.config:14-118`
- Schema: `nextflow_schema.json` (required for validation/help)
- Consume: in `workflows/deepumicaller.nf` or a subworkflow via `params.<name>`

**New Test:**
- Tests: `tests/main.nf.test` (add a `nextflow_pipeline` test block); fixtures in `tests/test_data/input/`, expected VCFs in `tests/test_data/expected_output/`

**New Config Profile:**
- Implementation: `conf/<profile>.config` + `includeConfig` entry in the `profiles` block of `nextflow.config:236-269`

## Special Directories

**`work/`:**
- Purpose: Nextflow task scratch directory
- Generated: Yes
- Committed: No

**`bin/__pycache__/`:**
- Purpose: Python bytecode cache
- Generated: Yes
- Committed: No (should be ignored)

**`tests/test_data/`:**
- Purpose: test fixtures and expected outputs
- Generated: No
- Committed: Yes

**`.planning/`:**
- Purpose: GSD planning documents (codebase maps, phase plans)
- Generated: By GSD tooling
- Committed: Yes

**`testing/`:**
- Purpose: currently empty; no files
- Generated: No
- Committed: Empty directory

---

*Structure analysis: 2026-10-10*
