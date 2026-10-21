---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# Coding Conventions

**Analysis Date:** 2026-10-10

## Languages in Scope

- **Nextflow (DSL2)** — pipeline logic: `main.nf`, `workflows/deepumicaller.nf`, `modules/local/**/main.nf`, `subworkflows/local/**/*.nf`
- **Python 3** — helper scripts executed inside processes: `bin/*.py`
- **Groovy** — test assertions inside `tests/main.nf.test`
- **YAML/JSON** — config and schema files: `conf/*.config`, `assets/schema_input.json`, `nextflow_schema.json`, `modules/local/*/meta.yml`

## Naming Patterns

**Files:**
- Python scripts: `snake_case.py` (e.g., `bin/check_samplesheet.py`, `bin/recompute_depth.py`, `bin/compute_umi_collisions.py`)
- One exception: `bin/VCF_Parser.py` (PascalCase, imported as a module by `bin/count_muts_per_cycle.py` via `from VCF_Parser import *`)
- Nextflow modules: `main.nf` inside a kebab-case or snake_case directory named after the tool/step (e.g., `modules/local/umicollisions/main.nf`, `modules/local/calling_vardict/calling_vardict_chunk/main.nf`)
- Config files: lowercase with `.config` extension in `conf/` (e.g., `conf/base.config`, `conf/modules.config`, `conf/results_outputs.config`)

**Nextflow processes:**
- `UPPER_SNAKE_CASE` process names (e.g., `UMICOLLISIONS`, `SAMPLESHEET_CHECK`, `CALLING_VARDICT_CHUNK`, `PATCH_DEPTH`, `FGUMI_CLIPBAM`)
- Local modules prefixed by tool family when wrapping a tool suite: `FGUMI_*` (e.g., `FGUMI_CLIPBAM` in `modules/local/fgumi/clipbam/main.nf`)
- Same module aliased at import site for reuse: `include { ALIGN_BAM as ALIGNRAWBAM } from '../modules/local/align_bam/main'` (`workflows/deepumicaller.nf:36-37`)

**Nextflow workflows:**
- `UPPER_SNAKE_CASE` (e.g., `DEEPUMICALLER` in `workflows/deepumicaller.nf`, `INPUT_CHECK` in `subworkflows/local/input_check.nf`, `POSTPROCESS_MUTATIONS` in `subworkflows/local/postprocess_mutations/main.nf`)

**Python functions:**
- `snake_case` (e.g., `check_samplesheet()`, `sniff_format()`, `parse_args()` in `bin/check_samplesheet.py`; `read_position_group_sizes()`, `compute_umi_collisions()` in `bin/compute_umi_collisions.py`)

**Python variables:**
- `snake_case`; module-level constants in `UPPER_SNAKE_CASE` (e.g., `requirementsDict` in `bin/check_samplesheet.py` is a legacy exception; `mutation_types`, `snv_color` in `bin/summarize_muts_per_cycle.py`; `refConvert`, `compBase` in `bin/count_muts_per_cycle.py`)

**Channels/variables in Nextflow:**
- `ch_` prefix for channels (e.g., `ch_multiqc_files`, `ch_ref_fasta` in `workflows/deepumicaller.nf`)
- `meta` map as first element of input/output tuples: `tuple val(meta), path(bam)`

## Code Style

**Formatting:**
- No enforced formatter for Python (no `pyproject.toml`, `setup.cfg`, `.flake8`, or `black`/`ruff` config in repo)
- No enforced formatter for Nextflow; style follows nf-core conventions visually (4-space indent, aligned `input:`/`output:` blocks)
- Markdown linting only: `.markdownlint.json` disables `MD013` (line length) and sets `MD024` to `siblings_only` (duplicate headings allowed across siblings)

**Linting:**
- No pre-commit hooks (`.pre-commit-config.yaml` absent)
- No nf-core lint CI (`.github/workflows/` absent — only `.github/ISSUE_TEMPLATE/` with `bug_report.md` and `feature_request.md`)
- No Python lint config detected; code style is by convention, not tooling

## Import Organization

**Python order (observed):**
1. Shebang `#!/usr/bin/env python` or `#!/usr/bin/env python3` (mixed usage; prefer `python3` for new scripts)
2. Module docstring
3. Stdlib imports (`os`, `sys`, `re`, `argparse`, `logging`, `collections`)
4. Third-party imports (`click`, `pandas`, `numpy`, `matplotlib`, `pysam`, `scipy`, `seaborn`)
5. Local imports (`from VCF_Parser import *` in `bin/count_muts_per_cycle.py`)

**Special pattern:** plotting scripts must force a non-interactive matplotlib backend before importing pyplot because they run headless in containers — see `bin/compute_umi_collisions.py:24-27`:

```python
import matplotlib
matplotlib.use("Agg")
```

**Nextflow imports:**
- Grouped under banner comment blocks: `IMPORT LOCAL MODULES/SUBWORKFLOWS`, then `IMPORT NF-CORE MODULES/SUBWORKFLOWS` (`workflows/deepumicaller.nf:11-95`)
- Aliasing at import for reuse of the same module: `include { FGUMI_CLIPBAM as CLIPBAM } ...` and `as CLIPBAMAM`
- Relative paths from the importing file (`../modules/local/...`, `../../modules/...` from subworkflows)

## CLI Parsing (Python) — Three Coexisting Styles

**Style 1: click (preferred for new scripts)** — used in `bin/recompute_depth.py`, `bin/compute_umi_collisions.py`, `bin/add_filter_nrich.py`, `bin/as_minus_xs.py`, `bin/family_size_plots.py`, `bin/summarize_muts_per_cycle.py`:

```python
@click.command()
@click.option('--mpileup-file', '-m', required=True, type=click.Path(exists=True), help='Path to the input mpileup file.')
@click.option('--output-filename', '-o', required=True, type=click.Path(), help='Path to the output VCF file.')
@click.option('--verbose', '-v', is_flag=True, default=False, help='Print verbose output.')
def cli(mpileup_file, vcf_file, output_filename, suffix):
    main(mpileup_file, vcf_file, output_filename, suffix)

if __name__ == '__main__':
    cli()
```

(`bin/recompute_depth.py:585-598`)

**Style 2: argparse (nf-core heritage scripts)** — `bin/check_samplesheet.py:270-303` uses `parse_args(argv=None)` with a `main(argv=None)` coordinator and `sys.exit(main())` entry point. `bin/count_muts_per_cycle.py:209` uses `ArgumentParser()` directly. `bin/filtervcf.py:157-203` and `bin/add_filter_from_bed.py:165-201` follow the same `parse_args(argv=None)` pattern with named flags (migrated from raw `sys.argv` on 2026-10-21; `filtervcf.py` uses a `--split-type` store_true flag, `add_filter_from_bed.py` uses `--vcf/--bed/--output/--filter-name`).

**Style 3: raw `sys.argv` (legacy, avoid in new code)** — `bin/compare_vcfs.py:39-43` is the only remaining user.

**Do this for new scripts:** click with `click.Path(exists=True)` for inputs, short+long flags, `help=` on every option, and a thin `cli()` wrapper delegating to a testable `main()` function.

## Error Handling

**Nextflow workflow level:**
- Fail fast with `log.error` + `exit 1` (`workflows/deepumicaller.nf:138-139`):

```groovy
} else {
    log.error "No reference FASTA was specified (--ref_fasta)."
    exit 1
}
```

- Input channel validation uses `exit 1, "message"` with the offending path echoed (`subworkflows/local/input_check.nf:47-56`):

```groovy
if (!file(row.fastq_1).exists()) {
    exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
}
```

**Nextflow process level:**
- Shell strict mode globally enabled (`nextflow.config:252`): `process.shell = ['/bin/bash', '-euo', 'pipefail']`
- Retry strategy on transient/kill exit codes (`conf/base.config:19-21`):

```groovy
errorStrategy = { task.exitStatus in [143,137,104,134,139,140] ? 'retry' : 'finish' }
maxRetries    = 3
```

- Memory-scaling on retry: `memory = { 4.GB * task.attempt }` throughout `conf/base.config`
- Selective `errorStrategy = 'ignore'` for non-critical QC/coverage processes (`conf/base.config:202`, `conf/modules.config:154,158`)
- In-script validation with explicit stderr + exit 1 (e.g., empty-VCF guard in `modules/local/calling_vardict/calling_vardict_chunk/main.nf:66-73`)

**Python level:**
- `sys.exit(1)` for validation failures, `sys.exit(2)` for missing input file (`bin/check_samplesheet.py:310`)
- `logger.critical(...)` followed by `sys.exit(1)` for row-level samplesheet errors (`bin/check_samplesheet.py:247-249`)
- `raise ValueError(...)` for structural problems (`bin/add_filter_from_bed.py:131-133`)
- `sys.stderr.write(...)` + `raise Exception` in third-party-derived code (`bin/count_muts_per_cycle.py:73-79`)
- Test utilities never throw: they skip malformed lines silently to keep comparisons robust (`bin/test_utilities/compare_vcfs_detailed.py:96-99`)

## Logging

**Framework:** mixed — no single standard.
- `logging` module with `[%(levelname)s] %(message)s` format and a `--log-level` flag: `bin/check_samplesheet.py:308`
- `click.echo()` for user-facing output in click-based scripts (e.g., `bin/recompute_depth.py`, `bin/as_minus_xs.py`)
- Plain `print()` in legacy scripts (`bin/add_filter_from_bed.py:163`)
- `sys.stderr.write()` for parser errors (`bin/count_muts_per_cycle.py:57-63`)

**Pattern for new scripts:** use `click.echo()` for progress, `sys.stderr` for errors; add a `--verbose`/`-v` flag when step-by-step output helps debugging.

## Comments

**When to Comment:**
- Banner blocks with `~~~~` or `====` rules delimit major sections in both Nextflow and Python files (e.g., `main.nf:1-8`, `tests/main.nf.test:1-27`)
- Inline comments explain *why* for non-obvious parameters (`nextflow.config:60-100` annotates every default)
- Provenance comments credit external origins: `bin/count_muts_per_cycle.py:2-5` and `bin/VCF_Parser.py:5-8` credit the Kennedy-Lab-UW Duplex-Seq-Pipeline
- `TODO`/`FIXME` markers used sparingly for known debt (e.g., `modules/local/fgumi/clipbam/main.nf:23`, `bin/recompute_depth.py:276`)

**Docstrings:**
- Google-style `Args:`/`Returns:` in nf-core-derived scripts (`bin/check_samplesheet.py:210-236`)
- NumPy-style `Parameters`/`Returns` sections with type annotations in pandas-heavy scripts (`bin/add_filter_from_bed.py:24-52`)
- Module-level docstrings describing inputs/outputs at the top of newer scripts (`bin/compute_umi_collisions.py:2-20`)
- Nextflow subworkflow inputs documented as trailing comments on `take:` blocks (`subworkflows/local/postprocess_mutations/main.nf:41-52`):

```groovy
take:
bam_n_index              // channel: [mandatory] [ val(meta), path (bam), pa
th (bamindex) ]
```

## Function Design

**Size:** No hard limit; largest script `bin/recompute_depth.py` (682 lines) keeps one `main()` plus focused helpers (`update_vcf_fields`, `remove_pstd_if_single`, `vartype`).

**Parameters:** keyword arguments with defaults for optional behavior (`verbose=False`, `output_file=None`, `suffix=''`).

**Return values:** pandas DataFrames for VCF-transforming helpers (`bin/add_filter_from_bed.py:52` returns `pd.DataFrame`); tuples for multi-value returns (`bin/count_muts_per_cycle.py:52` returns `(readlist, bqlist)`).

## Nextflow Process Conventions (prescriptive)

Every local module in `modules/local/*/main.nf` follows this skeleton — replicate it for new processes:

```nextflow
process UMICOLLISIONS {
    tag "$meta.id"                        // tag with sample id (or "$meta.id-${chunk_file.simpleName}" for chunked)
    label 'collect_duplex_metrics'        // resource label from conf/base.config

    container "docker.io/bbglab/deepcsa-core:0.0.1-alpha"   // or the conda+container pair below

    input:
    tuple val(meta), path(position_group_sizes), path(umi_counts), path(duplex_umi_counts)

    output:
    tuple val(meta), path("*.umi_collisions.tsv")   , emit: tsv
    path "versions.yml"                             , topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: ""
    prefix = "${meta.id}${prefix}"
    """
    compute_umi_collisions.py \\
        --sample-name ${meta.id} \\
        --output-file ${prefix}.umi_collisions.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.umi_collisions.tsv
    """
}
```

(`modules/local/umicollisions/main.nf`)

**Rules:**
1. **`tag`** — always `"$meta.id"`; add a suffix for chunked/sub-jobs (`modules/local/calling_vardict/calling_vardict_chunk/main.nf:2`)
2. **`label`** — always a resource label defined in `conf/base.config` (`alignment_intensive`, `bam_sort_compress`, `variant_calling`, `genomic_prep`, etc.). Never set cpus/memory inline in the module.
3. **Software provisioning** — either a single `container` directive (bbglab images like `docker.io/bbglab/deepcsa-core:0.0.1-alpha`) or the dual `conda` + singularity-aware `container` pair:

```nextflow
conda "bioconda::vardict-java=1.8.3"
container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0' :
    'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"
```

(`modules/local/calling_vardict/calling_vardict_chunk/main.nf:4-8`)
4. **`versions.yml`** — every process emits it with `topic: versions` and writes it via the `cat <<-END_VERSIONS` heredoc (26 modules do this)
5. **`stub:`** blocks — provide for expensive processes so nf-test/`stub` runs are cheap (9+ modules: `modules/local/patchdepth/main.nf:44-56`, `modules/local/expand_regions/main.nf:56-69`, etc.)
6. **`when:` guard** — optional `when: task.ext.when == null || task.ext.when` for skippable processes (`modules/local/calling_vardict/calling_vardict_chunk/main.nf:22-23`)
7. **`ext.args` pattern** — always `def args = task.ext.args ?: ''` and interpolate `${args}` so `conf/modules.config` can inject flags
8. **`ext.prefix` pattern** — always `def prefix = task.ext.prefix ?: ""` then `prefix = "${meta.id}${prefix}"` for output naming
9. **No `publishDir` in modules** — publishing is centralized in `conf/modules.config` and `conf/results_outputs.config` via `withName:` selectors (0 occurrences of `publishDir` in `modules/local/`). Default publish path derives from the process name: `"${params.outdir}/processing_files/${task.process.tokenize(':')[-1].tokenize('_')[0].toLowerCase()}"` (`conf/modules.config:12-16`)
10. **`meta.yml`** — nf-core modules carry `meta.yml`; most local modules do not (only some like `modules/local/align_bam/meta.yml`, `modules/local/count_ns/meta.yml`)

## Input Validation

- **Schema-based:** `assets/schema_input.json` (JSON Schema draft 2020-12) validates the samplesheet via the `nf-schema@2.3.0` plugin (`nextflow.config:283-285`); per-field `errorMessage` strings are mandatory
- **Pipeline params schema:** `nextflow_schema.json` with grouped `$defs` (nf-core layout)
- **Dual validation:** `bin/check_samplesheet.py` re-validates rows in Python (step-dependent required columns via `requirementsDict`, `bin/check_samplesheet.py:18-26`) because schema validation alone cannot express step-dependent requirements

## Config Conventions

- `nextflow.config` at root: defaults in `params {}` block with inline comments, `manifest` block, `plugins` block, profiles via `includeConfig` (`nextflow.config:236-243`)
- Profiles: `test`, `local`, `test_full`, `exome`, `irbcluster` — each maps to one file in `conf/`
- `conf/base.config`: resource labels in two parts — "PART 1: RESOURCE LABELS" and "PART 2: PROCESS-SPECIFIC OVERRIDES" with `withName: 'DEEPUMICALLER:PROCESS'` selectors
- `conf/modules.config`: per-module `ext.args`/`ext.prefix`/`publishDir` overrides
- `conf/results_outputs.config`: publish-path overrides only
- Institutional/reference data in `conf/igenomes.config` and `conf/general_files_IRB.config` (IRB cluster paths)

## Module Design

**Exports:** one process per `main.nf`; helper Groovy functions live in `modules/local/utils.nf` (`process_bams`, `clean_chr_names`) and are imported directly.

**Barrel files:** none; import paths point at specific `main.nf` files.

---

*Convention analysis: 2026-10-10*
