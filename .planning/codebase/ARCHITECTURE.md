---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
<!-- refreshed: 2026-10-10 -->

# Architecture

**Analysis Date:** 2026-10-10

## System Overview

```text
┌──────────────────────────────────────────────────────────────────────────────┐
│                          ENTRY POINT                                         │
│  `main.nf`  →  PIPELINE_INITIALISATION + DEEPUMICALLER                       │
└──────────────────────────────┬───────────────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                     ORCHESTRATION LAYER                                      │
│  `workflows/deepumicaller.nf`  (single 652-line workflow, params.step-gated) │
├──────────────────┬───────────────────┬──────────────────┬────────────────────┤
│  INPUT_CHECK     │  BAM_CALL_VARDICT │ POSTPROCESS_     │ VCF_ANNOTATE_ALL   │
│  (samplesheet)   │  _PARALLEL        │ MUTATIONS        │ (VEP annotation)   │
│  `subworkflows/  │  `subworkflows/   │ `subworkflows/   │ `subworkflows/     │
│   local/         │   local/bam_call_ │  local/          │  local/            │
│   input_check.nf`│   vardict_        │  postprocess_    │  vcf_annotate_all/ │
│                  │   parallel/`      │  mutations/`     │  main.nf`          │
└────────┬─────────┴─────────┬─────────┴────────┬─────────┴─────────┬──────────┘
         │                   │                  │                   │
         ▼                   ▼                  ▼                   ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│                       PROCESS LAYER (DSL2 modules)                           │
│  `modules/local/*`   custom processes (fgumi, VarDict, filters, metrics)     │
│  `modules/nf-core/*` vendored nf-core modules (fastp, fastqc, samtools, ...) │
│  `modules/local/utils.nf`  shared Groovy channel functions                   │
└──────────────────────────────┬───────────────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  TOOL LAYER                                                                  │
│  `bin/*` Python helpers (on PATH) · containers/conda per process             │
└──────────────────────────────┬───────────────────────────────────────────────┘
                               │
                               ▼
┌──────────────────────────────────────────────────────────────────────────────┐
│  OUTPUT (published via `conf/modules.config` + `conf/results_outputs.config`)│
│  ${params.outdir}/metrics · mutations_vcf · processing_files · pipeline_info │
└──────────────────────────────────────────────────────────────────────────────┘
```

## Component Responsibilities

| Component | Responsibility | File |
|-----------|----------------|------|
| `main.nf` | Entry point; wires initialisation + main workflow | `main.nf` |
| `PIPELINE_INITIALISATION` | Version banner, help text, param validation (nf-schema), config checks | `subworkflows/local/utils_nfcore_deepumicaller/main.nf` |
| `PIPELINE_COMPLETION` | onComplete hooks: email, summary, IM notification. **Invoked from `main.nf` (2026-10-21)** with email params + `DEEPUMICALLER.out.multiqc_report`; `hook_url` passed as `null` | `subworkflows/local/utils_nfcore_deepumicaller/main.nf` |
| `DEEPUMICALLER` | Main orchestration; all channel wiring and step gating | `workflows/deepumicaller.nf` |
| `INPUT_CHECK` | Samplesheet validation + meta-map construction per `params.step` | `subworkflows/local/input_check.nf` |
| `SAMPLESHEET_CHECK` | Schema validation of samplesheet CSV | `modules/local/samplesheet_check.nf` |
| `BAM_CALL_VARDICT_PARALLEL` | Chunked parallel VarDict calling (split → call → merge) | `subworkflows/local/bam_call_vardict_parallel/main.nf` |
| `POSTPROCESS_MUTATIONS` | Pileup, depth patching, region/N-rich/somatic filtering, muts-per-cycle | `subworkflows/local/postprocess_mutations/main.nf` |
| `VCF_ANNOTATE_ALL` | VEP annotation wrapper (snpEff/merge paths stubbed out) | `subworkflows/local/vcf_annotate_all/main.nf` |
| `process_bams`, `clean_chr_names` | Shared channel helpers for chromosome-split mode | `modules/local/utils.nf` |
| `bin/*.py` | Python utilities invoked inside processes (filters, plots, VCF parsing) | `bin/` |

## Pattern Overview

**Overall:** nf-core-style DSL2 pipeline — channel-based dataflow with meta-map tuples, resumable via `params.step` entry points.

**Key Characteristics:**
- Every process consumes/produces `tuple val(meta), path(...)` where `meta` is a Groovy map (`id`, `sample`, `parent_dna`, plus optional `read_structure`, `umi_file`, `chunk_name`).
- Modules are aliased at import time to reuse one process for multiple roles (e.g., `ALIGN_BAM as ALIGNRAWBAM` and `as ALIGNCONSENSUSBAM` in `workflows/deepumicaller.nf:36-37`).
- `params.step` gates `if` blocks so the pipeline can resume from: `mapping`, `preprocessed_crams`, `groupreadsbyumi`, `unmapped_consensus`, `allmoleculesfile`, `filterconsensus`, `calling`.
- Publishing is centralized in config (`conf/modules.config`, `conf/results_outputs.config`) via `withName:` blocks — modules contain no `publishDir` directives.
- Resource allocation via process `label` directives (`alignment_intensive`, `qc_memory_intensive`, `variant_calling`, ...) resolved in `conf/base.config` / `conf/exome.config`.
- Software versions collected through `topic: versions` outputs and `channel.topic('versions')` (`workflows/deepumicaller.nf:648`).

## Layers

**Entry/Init Layer:**
- Purpose: bootstrap, validate params, print help/logo
- Location: `main.nf`, `subworkflows/local/utils_nfcore_deepumicaller/main.nf`
- Depends on: nf-core utils subworkflows, `plugin/nf-schema`
- Used by: Nextflow runtime

**Orchestration Layer:**
- Purpose: define dataflow between processes; all channel topology lives here
- Location: `workflows/deepumicaller.nf`
- Contains: imports of all modules/subworkflows, `if (params.step ...)` gates, inline channel operators (`branch`, `join`, `groupTuple`, `combine`, `collectFile`)
- Depends on: all modules and local subworkflows
- Used by: nothing (terminal workflow)

**Subworkflow Layer (local):**
- Purpose: reusable multi-process groups with their own channel logic
- Location: `subworkflows/local/`
- Contains: `input_check.nf`, `bam_call_vardict_parallel/main.nf`, `postprocess_mutations/main.nf`, `vcf_annotate_all/main.nf`, `utils_nfcore_deepumicaller/main.nf`
- Depends on: `modules/local/*`, `modules/nf-core/*`, `subworkflows/nf-core/*`

**Process Layer:**
- Purpose: single-tool (or small pipeline) process definitions
- Location: `modules/local/<tool>/<process>/main.nf` (custom) and `modules/nf-core/<tool>/<process>/main.nf` (vendored)
- Contains: `process X { input/output/script }` with `conda` + `container` directives, `tag "$meta.id"`, `label`, `stub:` blocks
- Depends on: `bin/` scripts, containers/conda envs

**Tool/Script Layer:**
- Purpose: Python helpers automatically on PATH inside processes
- Location: `bin/` (e.g., `bin/filtervcf.py`, `bin/compute_umi_collisions.py`, `bin/as_minus_xs.py`, `bin/count_muts_per_cycle.py`)

## Data Flow

### Primary Request Path (params.step = 'mapping')

1. Samplesheet validated and split into meta-map channels (`subworkflows/local/input_check.nf:14-24`, meta built in `create_input_channel` at `subworkflows/local/input_check.nf:33`)
2. Targets BED expanded ±250bp (`modules/local/expand_regions/main.nf`, `EXPAND_PANEL` at `workflows/deepumicaller.nf:104-107`); converted to Picard interval list (`modules/nf-core/picard/bedtointervallist/main.nf`)
3. Optional adapter trim: `PRETRIMFASTQC` → `TRIMREADS` (fastp) → `FASTQC` (`workflows/deepumicaller.nf:139-158`)
4. Optional FASTQ splitting for parallelism: `SPLITFASTQ` + inline Groovy part-grouping (`workflows/deepumicaller.nf:161-186`)
5. FASTQ → unmapped BAM with UMI extraction: `FASTQTOBAM` (`modules/local/fgumi/fastqtobam/main.nf`); optional `CORRECTUMIS` (`modules/local/fgumi/correctumis/main.nf`)
6. Optional hard clipping: `TRIMBAM` (`modules/nf-core/bamutil/trimbam/main.nf`)
7. Alignment (fgumi fastq | bwa-mem3 | fgumi zipper | fgumi sort): `ALIGNRAWBAM` (`modules/local/align_bam/main.nf`)
8. Optional chromosome split: `SPLITBAMCHROM` + `process_bams()` helper (`modules/local/utils.nf:1-15`)
9. Template-coordinate sort: `SORTBAMRAWTEMPCOORDINATE` (`modules/local/fgumi/sort/main.nf`)
10. UMI grouping: `GROUPREADSBYUMI` (`modules/nf-core/fgumi/groupreadsbyumi/main.nf`)
11. Duplex metrics: `COLLECTSEQMETRICS` (+ on-target variant), `FAMILYMETRICS`, `UMICOLLISIONS` (`modules/local/familymetrics/main.nf`, `modules/local/umicollisions/main.nf`)
12. Consensus calling: `CALLCONSENSUSREADS` (`modules/nf-core/fgumi/callduplexconsensusreads/main.nf`)
13. Re-align consensus: `ALIGNCONSENSUSBAM` (same `modules/local/align_bam/main.nf`), optional `MERGEBAMCHROM` (`modules/local/mergebam/main.nf`)
14. All-molecules filtering: `ASMINUSXS` → `SAMTOOLSFILTERALLMOLECULES` → `SORTBAMAMFILTERED` (`modules/local/filter_reads/`)
15. Consensus filtering: parent_dna-grouped `MERGEBAMS` → `SORTBAMMERGED` → `FILTERCONSENSUSREADSAM`/`CLIPBAMAM`/`SORTBAMAMHQ` (all molecules) and `FILTERCONSENSUSREADSDUPLEX`/`CLIPBAM`/`SORTBAMDUPLEXCONS` (duplex) (`modules/local/fgumi/filterconsensusreads/main.nf`, `modules/local/fgumi/clipbam/main.nf`)
16. Depth + callable BED: `COMPUTEDEPTH` (`modules/nf-core/samtools/depth/main.nf`) → `CREATEBED` (`modules/local/createbed/main.nf`)
17. Variant calling: `BAM_CALL_VARDICT_PARALLEL` — `SPLIT_BED` → `CALLING_VARDICT_CHUNK` (parallel per chunk) → `MERGE_VARDICT_RESULTS` (`subworkflows/local/bam_call_vardict_parallel/main.nf:24-88`)
18. Post-processing: `POSTPROCESS_MUTATIONS` — mpileup (duplex + all-molecules), `NSXPOSITION`, `QUERYTABIX`, `PATCHDP`/`PATCHDPALL`, region filters (`FILTERLOWMAPPABLE`, `FILTERLOWCOMPLEX`, `FILTERNANOSEQSNP`, `FILTERNANOSEQNOISE`, `FILTERNRICH`), somatic filter (`FILTERVCFSOMATIC`), `MUTSPERPOS` + `COHORTMUTSPERPOS` (`subworkflows/local/postprocess_mutations/main.nf:41-227`)
19. Optional annotation: `VCFANNOTATE` → `VCF_ANNOTATE_ENSEMBLVEP` → `ENSEMBLVEP_VEP` (`subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf`)
20. Mutation signatures: `SIGPROFPLOT`/`SIGPROFPLOTPUR`/`SIGPROFPLOTPYR` (`modules/local/sigprofiler/matrixgenerator/main.nf`)
21. Reporting: `CUSTOM_DUMPSOFTWAREVERSIONS` + `MULTIQC`/`MULTIQCDUPLEX` (`workflows/deepumicaller.nf:648-676`)

### Resume/Entry-Point Flow (params.step ≠ 'mapping')

Each `if (params.step in [...])` block reassigns its input channel from `INPUT_CHECK.out.reads` when the step is the entry point:
- `preprocessed_crams` → `RETAGFROMCRAM` (`modules/local/fgumi/retag/main.nf`)
- `groupreadsbyumi` → `pre_consensus_bams = INPUT_CHECK.out.reads` (`workflows/deepumicaller.nf:322-325`)
- `unmapped_consensus` → `UNMAPBAM` (`modules/local/unmap/main.nf`)
- `allmoleculesfile` → `all_molecules_bam_complete_n_index = INPUT_CHECK.out.reads` (`workflows/deepumicaller.nf:432-434`)
- `filterconsensus` → `duplex_filtered_init_bam = INPUT_CHECK.out.reads` (`workflows/deepumicaller.nf:461-463`)
- `calling` → `cons_duplex_bam` / `bam_n_index_duplex_clean` from `INPUT_CHECK.out.reads` (`workflows/deepumicaller.nf:565-569`)

Samplesheet columns differ per step (`fastq_1/fastq_2/read_structure`, `preprocessed_cram`, `duplexbam/csi`, `bam`) — see `create_input_channel` in `subworkflows/local/input_check.nf:33-84`.

### Parallel VarDict Flow

1. `SPLIT_BED` chunks the callable BED (`params.vardict_chunks`, default 2) (`modules/local/calling_vardict/split_bed/main.nf`)
2. Chunks transposed/combined with BAM+BAI by `meta.id`; `chunk_name` added to meta (`subworkflows/local/bam_call_vardict_parallel/main.nf:36-50`)
3. `CALLING_VARDICT_CHUNK` runs VarDict per chunk (`modules/local/calling_vardict/calling_vardict_chunk/main.nf`)
4. Chunk VCFs/TSVs grouped back by clean meta and merged (`MERGE_VARDICT_RESULTS`, `modules/local/calling_vardict/merge_vardict_results/main.nf`)

**State Management:**
- No persistent state; all state is in channel topology and the `meta` maps.
- Aggregated summaries written directly with `collectFile(..., storeDir: "${params.outdir}/...")` (e.g., `workflows/deepumicaller.nf:352-357`, `:447-451`).
- Split-by-chromosome mode reconstructs sample-level files via `clean_chr_names()` (`modules/local/utils.nf:17-32`).

## Key Abstractions

**Meta Map:**
- Purpose: per-sample/per-item metadata threaded through every channel
- Examples: built in `subworkflows/local/input_check.nf:33-84`; cloned/extended throughout `workflows/deepumicaller.nf`
- Pattern: Groovy LinkedHashMap; `meta.id` unique per item, `meta.sample` groups lanes, `meta.parent_dna` groups biological replicates for merging

**Module Aliasing:**
- Purpose: reuse one process definition for multiple pipeline roles
- Examples: `workflows/deepumicaller.nf:36-37` (`ALIGN_BAM` ×2), `:44-49` (`FGUMI_CLIPBAM`, `FGUMI_FILTERCONSENSUSREADS`, `FGUMI_COLLECTDUPLEXSEQMETRICS`, `SIGPROFILER_MATRIXGENERATOR` ×3 each)
- Pattern: `include { PROC as ALIAS } from '...'`; per-alias behavior tuned via `task.ext.args` in `conf/modules.config`

**Topic Channels (versions):**
- Purpose: collect `versions.yml` from every process without explicit wiring
- Examples: `path "versions.yml", topic: versions` in all modules; consumed at `workflows/deepumicaller.nf:648`
- Pattern: nf-core `topic:` output + `channel.topic('versions')`

**Shared Channel Functions:**
- Purpose: chromosome-split bookkeeping
- Examples: `process_bams()`, `clean_chr_names()` in `modules/local/utils.nf`
- Pattern: plain Groovy functions imported into the workflow

## Entry Points

**`main.nf`:**
- Location: `main.nf`
- Triggers: `nextflow run .` (manifest `mainScript` in `nextflow.config:247`)
- Responsibilities: run `PIPELINE_INITIALISATION(params.version, params.help, params.validate_params, params.monochrome_logs, args, params.outdir)` then `DEEPUMICALLER()`

**`tests/main.nf.test`:**
- Location: `tests/main.nf.test`
- Triggers: `nf-test test` (config `nf-test.config`, `tests/nextflow.config`)
- Responsibilities: pipeline-level functional tests against `tests/test_data/input/*.csv` with expected VCFs in `tests/test_data/expected_output/`

## Architectural Constraints

- **DSL:** `nextflow.enable.dsl = 2` (`main.nf:8`); Nextflow `!>=25.04.1` (`nextflow.config:249`); plugin `nf-schema@2.3.0` (`nextflow.config:283`)
- **Single-workflow orchestration:** all topology lives in `workflows/deepumicaller.nf` (652 lines); adding a step means editing this file's `if`-gated blocks.
- **Meta-map contract:** downstream `join`/`groupTuple` operations depend on exact `meta.id`/`meta.sample`/`meta.parent_dna` semantics; changing key names breaks joins (e.g., chromosome-split suffixing `${meta.id}_${chrom}` in `modules/local/utils.nf:12`).
- **Reference requirement:** `params.ref_fasta` mandatory with `.dict` and `.amb` indexes co-located (`workflows/deepumicaller.nf:96-106`); alignment locates the FASTA via `find ... *.amb` (`modules/local/align_bam/main.nf:33-34`).
- **Config load order:** `nextflow.config` → `conf/base.config` → profiles → `conf/modules.config` → `conf/results_outputs.config` (`nextflow.config:150, 271, 311-312`).
- **Containers/conda pinned per module:** each module declares its own `container`/`conda`; several custom images under `docker.io/ferriolcalvet/*` and `docker.io/bbglab/*`.
- **Publishing centralized:** `publishDir` only in `conf/modules.config` and `conf/results_outputs.config`; `versions.yml` filtered via `saveAs`.

## Anti-Patterns

### Monolithic Workflow File

**What happens:** All channel wiring, branching, and step gating for the entire pipeline is in one 652-line workflow (`workflows/deepumicaller.nf`).
**Why it's wrong:** Hard to navigate and test in isolation; inline Groovy (e.g., splitfastq part grouping at `workflows/deepumicaller.nf:163-186`, chromosome regrouping at `:415-437`) duplicates logic that belongs in subworkflows/functions.
**Do this instead:** Extract each `if (params.step ...)` block into a local subworkflow under `subworkflows/local/` (pattern already used by `bam_call_vardict_parallel`).

### Commented-Out Code Paths

**What happens:** snpEff/merge annotation branches are fully commented in `subworkflows/local/vcf_annotate_all/main.nf:17-27`; pysam container alternatives commented in `modules/local/mutsperpos/compute/main.nf:4-10`.
**Why it's wrong:** Dead code obscures the actual annotation path (VEP only).
**Do this instead:** Remove or gate behind params like the existing `params.annotate_mutations` flag.

### Hardcoded Cluster Paths in Profiles

**What happens:** `conf/test.config` and `conf/general_files_IRB.config` reference `/data/bbg/datasets/...` absolute paths.
**Why it's wrong:** Tests and the IRB profile fail off-cluster; test data is not self-contained.
**Do this instead:** Ship small test fixtures under `tests/test_data/input/` (already exists) and point `conf/test.config` at them.

## Error Handling

**Strategy:** Process-level retry on transient signals; fail-fast on invalid input.

**Patterns:**
- `errorStrategy = { task.exitStatus in [143,137,104,134,139,140] ? 'retry' : 'finish' }` with `maxRetries 3` (`conf/base.config:20-22`); memory scaled by `task.attempt`.
- Samplesheet file existence checked with `exit 1, "ERROR: ..."` per column (`subworkflows/local/input_check.nf:44-56`).
- Missing reference: `log.error "No reference FASTA was specified (--ref_fasta)."; exit 1` (`workflows/deepumicaller.nf:105-107`).
- Optional QC/filters bypassed via channel pass-through when params unset (e.g., `ch_after_lowmap = ch_vcf_with_bed` in `subworkflows/local/postprocess_mutations/main.nf:110-112`).
- Temporary process ignoring via `conf/tmp_quick_fixes.config` (`errorStrategy = 'ignore'` for `.*COVERAGE.*`).

## Cross-Cutting Concerns

**Logging:** Nextflow native logging; banner/logo via `logo()` in `subworkflows/local/utils_nfcore_deepumicaller/main.nf:57`; completion email/summary via `PIPELINE_COMPLETION` hooks, invoked from `main.nf` (2026-10-21). IM/webhook notification exists upstream but `hook_url` is passed as `null` (unwired).

**Validation:** Samplesheet validated against `assets/schema_input.json` by `SAMPLESHEET_CHECK` (`modules/local/samplesheet_check.nf`); pipeline params validated by `nf-schema` against `nextflow_schema.json` (`subworkflows/local/utils_nfcore_deepumicaller/main.nf:63-71`).

**Authentication:** None — no external service auth; references resolved from local paths or `params.igenomes_base` (`conf/igenomes.config`).

**Reporting:** MultiQC twice (mid-run `MULTIQCDUPLEX` when `params.perform_qcs`, final `MULTIQC`); execution timeline/report/trace/DAG auto-enabled into `${params.outdir}/pipeline_info/` (`nextflow.config:262-296`).

---

*Architecture analysis: 2026-10-10*
