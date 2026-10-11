---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# Technology Stack

**Analysis Date:** 2026-10-10

## Languages

**Primary:**
- **Nextflow (DSL2)** - Pipeline orchestration language. Entry point `main.nf` (`nextflow.enable.dsl = 2`). All workflows in `workflows/deepumicaller.nf`, subworkflows in `subworkflows/local/` and `subworkflows/nf-core/`, processes in `modules/local/` and `modules/nf-core/`.
- **Python 3** - Helper scripts in `bin/` (e.g., `bin/check_samplesheet.py`, `bin/count_muts_per_cycle.py`, `bin/compute_umi_collisions.py`, `bin/family_size_plots.py`, `bin/as_minus_xs.py`, `bin/filtervcf.py`, `bin/VCF_Parser.py`). Scripts use shebangs `#!/usr/bin/env python` or `#!/usr/bin/env python3` and are executed inside containers, not on the host.
- **Groovy** - Nextflow config files (`nextflow.config`, `conf/*.config`) and the sendmail template rendering in `subworkflows/nf-core/utils_nfcore_pipeline/main.nf`.

**Secondary:**
- **Shell/Bash** - Process scripts inside `.nf` modules; `process.shell = ['/bin/bash', '-euo', 'pipefail']` enforced in `nextflow.config`.

## Runtime

**Environment:**
- **Nextflow >= 25.04.1** (enforced via `nextflowVersion = '!>=25.04.1'` in the `manifest` block of `nextflow.config`)
- DSL2 (`nextflow.enable.dsl = 2` in `main.nf`)
- Java (bundled with Nextflow) for VarDict-Java and Picard execution inside containers

**Package Manager:**
- **Conda/Mamba** - Optional via `conda` and `mamba` profiles in `nextflow.config` (per-process `conda` directives in modules)
- **Containers (primary)** - Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer profiles all defined in `nextflow.config`; Singularity is the profile used for testing (`nf-test.config` sets `profile "test,singularity"`)
- Lockfile: Not applicable (no conda `environment.yml` lock; dependencies pinned per-process via `conda`/`container` directives)

**Nextflow Plugins:**
- `nf-schema@2.3.0` - Parameter validation and samplesheet channel creation (declared in `plugins` block of `nextflow.config`; used via `include { paramsSummaryMap } from 'plugin/nf-schema'` in `workflows/deepumicaller.nf`)
- `nft-utils@0.0.3` - Loaded in `nf-test.config` for pipeline-level tests

## Frameworks

**Core:**
- **nf-core modules framework** - 14 pinned nf-core modules under `modules/nf-core/` (see `modules.json` for git SHAs from `https://github.com/nf-core/modules.git`, branch `master`)
- **Local DSL2 modules** - ~30 process modules under `modules/local/` following the nf-core module structure (`main.nf` + `meta.yml`)
- Pipeline is a fork of **nf-core/fastquorum** (per `CHANGELOG.md`: "Initial fork of nf-core/fastquorum, created with the nf-core template"), rebranded as `bbglab/deepUMIcaller` v1.1.0 "Montseny"

**Testing:**
- **nf-test** - Pipeline-level tests in `tests/main.nf.test`; config in `nf-test.config` (testsDir `tests`, workDir from `DEEPUMICALLER_TEST_WORKDIR` env or `.nf-test/`, ignores `modules/nf-core/**` and `subworkflows/nf-core/**`)
- Test profiles: `test` (`conf/test.config`), `test_full` (`conf/test_full.config`)

**Build/Dev:**
- No compiled build step; the pipeline is pure Nextflow + containerized tools
- `markdownlint` config present (`.markdownlint.json`)

## Key Dependencies

**Critical (bioinformatics tools, pinned via container/conda directives):**

| Tool | Version | Used for | Module |
|------|---------|----------|--------|
| fgumi (fgbio fork) | 0.7.0 (conda `bioconda::fgumi=0.7.0`); unversioned `bioconda::fgumi` in `modules/local/fgumi/*` | UMI grouping, consensus calling, duplex metrics | `modules/local/fgumi/*`, `modules/local/align_bam/main.nf` |
| bwa-mem3 | 0.14.0 (`bioconda::bwa-mem3=0.14.0`) | Alignment | `modules/local/align_bam/main.nf` |
| VarDict-Java | 1.8.3 (`bioconda::vardict-java=1.8.3`) | Variant calling | `modules/local/calling_vardict/*` |
| samtools | 1.20 (local modules), 1.16.1/1.17/1.24 (nf-core modules) | BAM manipulation, mpileup, depth | `modules/local/mergebam/main.nf`, `modules/local/splitbamchrom/main.nf`, `modules/local/unmap/main.nf`, `modules/nf-core/samtools/*` |
| Ensembl VEP | 111.0 (also 102.0/108.2 variants present in `modules/nf-core/ensemblvep/`) | Variant annotation (offline cache mode) | `modules/nf-core/ensemblvep/vep/main.nf` |
| SigProfilerMatrixGenerator | 1.3.5 (container `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5`) | Mutational signature matrices/plots | `modules/local/sigprofiler/matrixgenerator/main.nf` |
| fastp | 0.23.2 | Read trimming | `modules/nf-core/fastp/main.nf` |
| FastQC | 0.12.1 | Read QC | `modules/nf-core/fastqc/main.nf` |
| Picard | 3.0.0 / 3.3.0 | BedToIntervalList, etc. | `modules/nf-core/picard/*` |
| Qualimap | 2.2.2d | BAM QC | `modules/nf-core/qualimap/bamqc/main.nf` |
| bamutil | 1.0.15 | trimBam (clipping) | `modules/nf-core/bamutil/trimbam/main.nf` |
| bedtools | 2.31.0 / 2.31.1 | BED operations | `modules/local/bedtools/merge/main.nf`, `modules/nf-core/bedtools/*` |
| tabix | 1.11 | Indexed VCF/BED queries | `modules/local/count_ns/main.nf`, `modules/local/tabix_mpileup/main.nf`, `modules/local/expand_regions/main.nf` |
| seqkit | 2.3.1 | FASTQ splitting | `modules/local/splitfastq/main.nf` |
| pybedtools | 0.9.1 | BED-based VCF filtering | `modules/local/filter/from_bed/main.nf` |
| SnpEff | (module present: `modules/nf-core/snpeff/` per `modules.json`; annotation path currently commented out in `subworkflows/local/vcf_annotate_all/main.nf`) | Alternative annotation | `modules/nf-core/snpeff/*` |

**In-house / bbglab containers (not on Bioconda):**

| Container | Used by |
|-----------|---------|
| `docker.io/bbglab/deepcsa-core:0.0.1-alpha` | `modules/local/familymetrics/main.nf`, `modules/local/patchdepth/main.nf`, `modules/local/umicollisions/main.nf` |
| `docker.io/bbglab/pysam-0.23.3:latest` | `modules/local/filter_reads/asminusxs/main.nf`, `modules/local/filter/nrich/main.nf` |
| `docker.io/ferriolcalvet/pysam` | `modules/local/mutsperpos/compute/main.nf` |
| `docker.io/ferriolcalvet/fgumi:v-retag` | `modules/local/fgumi/retag/main.nf` |
| `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5` | `modules/local/sigprofiler/matrixgenerator/main.nf` |
| `community.wave.seqera.io/library/bwa-mem3_fgumi:ab0c353b848ba934` | `modules/local/align_bam/main.nf` (Wave-built fused container) |
| `biocontainers/seaborn:0.12.2_cv1` | `modules/local/filtervcf/main.nf` |

**Python libraries used by `bin/` scripts (provided by the containers above):**
- `pysam` (0.21.0 / 0.23.3 depending on module)
- `pandas`, `numpy`, `matplotlib`, `seaborn`, `scipy`, `click`

## Configuration

**Environment:**
- All runtime configuration via Nextflow `params` (defaults in `nextflow.config`), validated by `nf-schema@2.3.0` against `nextflow_schema.json` (param groups in `$defs`: `input_output_options`, `main_options`, `raw_data_processing_options`, `mutation_filtering_options`, `asxs_filtering_params`, `grouping`, `consensus_calling`, `consensus_filtering_general`, `all_molecules_filtering`, `calling_parameters`, `consensus_filtering`, `ensembl_vep_options`, `read_processing_options`, `output_naming_options`, `institutional_config_options`, `max_job_request_options`, `generic_options`)
- Input samplesheet validated against `assets/schema_input.json` (fields: `sample`, `fastq_1`, `fastq_2`, `umi_file`, `bam`, `duplexbam`, `preprocessed_cram`, `csi`, `read_structure`, `parent_dna`)
- No `.env` files; secrets not used (pipeline runs on HPC with local references)

**Build/Config files:**
- `nextflow.config` - Global defaults, profiles, manifest, plugins, timeline/report/trace/DAG reporting
- `conf/base.config` - Default process resources and retry strategy (`errorStrategy` retry on exit codes 143/137/104/134/139/140, `maxRetries 3`)
- `conf/modules.config` - Per-process `ext.args`, `ext.prefix`, `publishDir` overrides
- `conf/results_outputs.config` - Publish directory routing
- `conf/igenomes.config` - iGenomes reference definitions (`params.igenomes_base` base path; GRCh37, GRCh38, GRCm38, TAIR10, EB2, UMD3.1, WBcel235, CanFam3.1, and more)
- `conf/local.config` - Local executor resource limits (16 CPUs, 60 GB)
- `conf/exome.config` - Exome-optimized resource labels
- `conf/general_files_IRB.config` - IRB cluster site config (VEP cache, reference FASTA, masking BED files, Singularity cache dirs)
- `conf/test.config`, `conf/test_full.config` - Test profiles
- `conf/tmp_quick_fixes.config` - Ad-hoc overrides
- `nf-test.config` - nf-test runner config
- `tower.yml` - Seqera Platform (Tower) report file display config

**Profiles defined in `nextflow.config`:**
`debug`, `conda`, `mamba`, `docker` (registry `quay.io`), `arm`, `singularity`, `podman`, `shifter`, `charliecloud`, `apptainer`, `gitpod`, `test`, `local`, `test_full`, `exome`, `irbcluster`

## Platform Requirements

**Development:**
- Nextflow >= 25.04.1 with a container engine (Docker or Singularity/Apptainer)
- nf-test binary for running `tests/main.nf.test`
- For IRB cluster: Singularity cache at `/data/bbg/datasets/pipelines/nextflow_containers` (set in `conf/general_files_IRB.config`)

**Production:**
- HPC cluster (SLURM-style) or local execution; `conf/base.config` assumes software available via containers
- Reference data: iGenomes base (`--igenomes_base`) or custom FASTA via `--ref_fasta`; VEP cache via `--vep_cache` (default at IRB: `/data/bbg/datasets/vep`)
- Resource ceilings: `max_memory 1024.GB`, `max_cpus 192`, `max_time 240.h` (`nextflow.config`)

---

*Stack analysis: 2026-10-10*
