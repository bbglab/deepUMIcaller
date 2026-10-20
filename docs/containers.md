# Containers used

This document inventories every container used by the pipeline, where it is
defined, and — when known — how it was generated.

Containers are defined in three places:

1. **Per-process `container` directives** in `modules/local/**/main.nf` and
   `modules/nf-core/**/main.nf`.
2. **Label-based overrides** in `conf/modules.config` (`withLabel:` blocks).
3. Most `nf-core`-style modules use the conditional pattern
   `workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container`
   to pick a Galaxy Project pre-built Singularity image
   (`https://depot.galaxyproject.org/singularity/...`) instead of pulling the
   Docker image from BioContainers.

---

## Custom containers

### Wave (Seqera) built containers

Containers built on demand with [Wave](https://wave.seqera.io/), which assembles
a container from the conda dependencies listed in the module (`conda`
directive) and hosts it on `community.wave.seqera.io`.

| Container | Used by | Conda spec | Build page |
|---|---|---|---|
| `community.wave.seqera.io/library/bwa-mem3_fgumi:6069fbf2f2bcec9f` | `ALIGN_BAM` (`modules/local/align_bam/main.nf`) | `bioconda::fgumi=0.7.0 bioconda::bwa-mem3=0.14.0` | https://wave.seqera.io/view/builds/bd-6069fbf2f2bcec9f_1 |
| community.wave.seqera.io/library/fgumi_samtools:f4eabe7668db9d70 | `retag` | `bioconda::fgumi=0.8.0 bioconda::samtools=1.24` | https://wave.seqera.io/view/builds/bd-f4eabe7668db9d70_1 |

### BBGLab Docker Hub containers (`docker.io/bbglab/...`)

Built in-house by BBGLab. Recipes live in the
[bbglab/containers-recipes](https://github.com/bbglab/containers-recipes)
repository and images are published on
[Docker Hub `bbglab`](https://hub.docker.com/u/bbglab).

| Container | Used by | Origin / recipe |
|---|---|---|
| `docker.io/bbglab/deepcsa-core:0.0.1-alpha` | `FAMILYSIZEMETRICS` (`modules/local/familymetrics`), `UMICOLLISIONS` (`modules/local/umicollisions`), `PATCH_DEPTH` (`modules/local/patchdepth`), `SUMMARIZE_MUTS_PER_POS` (`modules/local/mutsperpos/summarize`) | Recipe [`deepcsa_core/`](https://github.com/bbglab/containers-recipes/tree/main/deepcsa_core) in `bbglab/containers-recipes`. Built from `mambaorg/micromamba:2.0.5-debian12-slim` with an explicit conda environment (`deepcsa.txt`) plus `requirements.txt` (includes `pysam==0.23.0`, `pandas`, `matplotlib`, `seaborn`, ...), and pre-downloaded `bgdata` genome references (hg38, mm39). See its [README](https://github.com/bbglab/containers-recipes/blob/main/deepcsa_core/README.md) for the rebuild workflow. |
| `docker.io/bbglab/pysam-0.23.3:latest` | `ASMINUSXS` (`modules/local/filter_reads/asminusxs`), `FILTER_N_RICH` (`modules/local/filter/nrich`) | No dedicated recipe found in `bbglab/containers-recipes`. The closest recipe is [`mutsperpos/Dockerfile`](https://github.com/bbglab/containers-recipes/tree/main/mutsperpos) (Ubuntu 22.04 + `pysam==0.23.3` + numpy/pandas/matplotlib), which may be the same or a sibling build. **Origin to confirm.** |

### Personal Docker Hub containers (`docker.io/ferriolcalvet/...`)

Built and pushed by the pipeline maintainer to a personal Docker Hub account.

| Container | Used by | Origin / recipe |
|---|---|---|
| `docker.io/ferriolcalvet/pysam` | `MUTS_PER_POS` (`modules/local/mutsperpos/compute`) | Likely built from the [`mutsperpos/Dockerfile`](https://github.com/bbglab/containers-recipes/tree/main/mutsperpos) recipe in `bbglab/containers-recipes` (Ubuntu 22.04, `pysam==0.23.3`, numpy, pandas, matplotlib), but the image is hosted under the personal account and the exact build is not traceable. **Origin to confirm.** |
| `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5` | `SIGPROFILER_MATRIXGENERATOR` (`modules/local/sigprofiler/matrixgenerator`) | Closest recipe: [`sigprofiler_matrix/Dockerfile`](https://github.com/bbglab/containers-recipes/tree/main/sigprofiler_matrix) in `bbglab/containers-recipes` (`python:3.9.16-slim-bullseye` + `SigProfilerMatrixGenerator` + `click`, with GRCh37/GRCh38/mm39/mm10 genomes pre-installed). The pinned `1.3.5` personal build is not documented in the recipes repo. **Origin to confirm.** |
| `docker.io/ferriolcalvet/fgumi:v-retag` | `FGUMI_RETAGFROMCRAM` (`modules/local/fgumi/retag`) | No recipe found in `bbglab/containers-recipes` (no `fgumi` recipe exists there). Presumably a custom image bundling `fgumi` + `samtools=1.24` (matching the module's `conda "bioconda::fgumi bioconda::samtools=1.24"`). **Origin to confirm.** |

---

## Label-based containers (`conf/modules.config`)

These are applied to every process carrying the corresponding `label`, which
overrides any `container` directive in the module itself.

| Label | Container | Used by |
|---|---|---|
| `fgumi_tools` | `quay.io/biocontainers/fgumi:0.7.0--h54198d6_0` | `FGUMI_CLIPBAM`, `FGUMI_COLLECTDUPLEXSEQMETRICS`, `FGUMI_CORRECTUMIS`, `FGUMI_FASTQTOBAM`, `FGUMI_FILTERCONSENSUSREADS`, `FGUMI_RETAGFROMCRAM`*, `FGUMI_SORT` (all `modules/local/fgumi/*`), `FGUMI_GROUPREADSBYUMI`, `FGUMI_CALLDUPLEXCONSENSUSREADS` (`modules/nf-core/fgumi/*`) |
| `samtools_module` | `https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0` | **None currently** — no process in the codebase carries the `samtools_module` label. This override is dead configuration and can be removed or re-purposed. |

\* `FGUMI_RETAGFROMCRAM` also declares its own `container` directive
(`docker.io/ferriolcalvet/fgumi:v-retag`); the label override takes precedence,
so the per-module container is effectively unused.

The `fgumi` biocontainers image is the standard BioContainers build of the
bioconda package (`fgumi=0.7.0`), matching the `conda "bioconda::fgumi"`
directive present in all these modules.

---

## Standard BioContainers / Galaxy images

Standard community images built from bioconda recipes, published on
[BioContainers](https://biocontainers.pro/) (quay.io/biocontainers, pulled via
Docker Hub `biocontainers/...`) with pre-built Singularity mirrors at
[depot.galaxyproject.org](https://depot.galaxyproject.org/singularity). These
require no custom build; the version tags map 1:1 to bioconda package builds.

| Container (Docker / Singularity) | Used by | Module |
|---|---|---|
| `biocontainers/python:3.8.3` / `depot.galaxyproject.org/singularity/python:3.8.3` | `SAMPLESHEET_CHECK` | `modules/local/samplesheet_check.nf` |
| `biocontainers/bedtools:2.31.0--hf5e1c6e_2` / `.../bedtools:2.31.0--hf5e1c6e_2` | `BEDTOOLS_MERGE`, `CREATEBED` | `modules/local/bedtools/merge`, `modules/local/createbed` |
| `biocontainers/bedtools:2.31.1--hf5e1c6e_0` / `.../bedtools:2.31.1--hf5e1c6e_0` | `BEDTOOLS_COVERAGE` | `modules/nf-core/bedtools/coverage` |
| `biocontainers/tabix:1.11--hdfd78af_0` / `.../tabix:1.11--hdfd78af_0` | `COUNT_NS`, `EXPAND_REGIONS`, `QUERY_TABIX` | `modules/local/count_ns`, `modules/local/expand_regions`, `modules/local/tabix_mpileup` |
| `biocontainers/samtools:1.20--h50ea8bc_1` / `.../samtools:1.20--h50ea8bc_1` | `SAMTOOLS_FILTER`, `MERGEBAM`, `SPLITBAMCHROM`, `UNMAP`, `SAMTOOLS_DEPTH` | `modules/local/filter_reads/samtools`, `modules/local/mergebam`, `modules/local/splitbamchrom`, `modules/local/unmap`, `modules/nf-core/samtools/depth` |
| `biocontainers/samtools:1.17--h00cdaf9_0` / `.../samtools:1.17--h00cdaf9_0` | `SAMTOOLS_MPILEUP` | `modules/nf-core/samtools/mpileup` |
| `biocontainers/vardict-java:1.8.3--hdfd78af_0` / `.../vardict-java:1.8.3--hdfd78af_0` | `CALLING_VARDICT_CHUNK`, `MERGE_VARDICT_RESULTS` | `modules/local/calling_vardict/*` |
| `biocontainers/pybedtools:0.9.1--py38he0f268d_0` / `.../pybedtools:0.9.1--py38he0f268d_0` | `FILTER_FROM_BED` | `modules/local/filter/from_bed` |
| `biocontainers/seqkit:2.3.1--h9ee0642_0` / `.../seqkit:2.3.1--h9ee0642_0` | `SPLITFASTQ` | `modules/local/splitfastq` |
| `biocontainers/fastp:0.23.2--h79da9fb_0` / `.../fastp:0.23.2--h79da9fb_0` | `FASTP` | `modules/nf-core/fastp` |
| `biocontainers/fastqc:0.12.1--hdfd78af_0` / `.../fastqc:0.12.1--hdfd78af_0` | `FASTQC` | `modules/nf-core/fastqc` |
| `biocontainers/multiqc:1.29--pyhdfd78af_0` / `.../multiqc:1.29--pyhdfd78af_0` | `MULTIQC` | `modules/nf-core/multiqc` |
| `biocontainers/multiqc:1.18--pyhdfd78af_0` / `.../multiqc:1.18--pyhdfd78af_0` | `DUMPSOFTWAREVERSIONS` (reuses the MultiQC image because it needs `pyyaml`) | `modules/nf-core/custom/dumpsoftwareversions` |
| `biocontainers/picard:3.0.0--hdfd78af_1` / `.../picard:3.0.0--hdfd78af_1` | `PICARD_BEDTOINTERVALLIST` | `modules/nf-core/picard/bedtointervallist` |
| `biocontainers/picard:3.3.0--hdfd78af_0` / `.../picard:3.3.0--hdfd78af_0` | `PICARD_MERGESAMFILES` | `modules/nf-core/picard/mergesamfiles` |
| `biocontainers/qualimap:2.2.2d--1` / `.../qualimap:2.2.2d--1` | `QUALIMAP_BAMQC` | `modules/nf-core/qualimap/bamqc` |
| `biocontainers/bamutil:1.0.15--h2e03b76_1` / `.../bamutil:1.0.15--h2e03b76_1` | `TRIMBAM` | `modules/nf-core/bamutil/trimbam` |
| `biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0` (default, `vep_cache_version = 111`) or `108.2--pl5321h4a94de4_0` / `102.0--pl526hecda079_0` depending on `params.vep_cache_version` | `ENSEMBLVEP_VEP` | `modules/nf-core/ensemblvep/vep` |

### Semi-custom BioContainers image

| Container | Used by | Notes |
|---|---|---|
| `biocontainers/seaborn:0.12.2_cv1` | `FILTERMUTATIONS` (`modules/local/filtervcf`) | The `_cv1` suffix indicates a **custom variant build**, not a standard bioconda-derived BioContainers image. There is no Galaxy Singularity mirror for it (the module carries a TODO: *"update this in the nfcore format once the container is available in biocontainers and galaxy singularity"*). A related recipe exists in `bbglab/containers-recipes` under [`filtermuts/`](https://github.com/bbglab/containers-recipes/tree/main/filtermuts) (a `Dockerfile2222` file declares `seaborn 0.12.2` labels), suggesting this image was built for the filtermuts tooling. **Origin to confirm.** |

---

## Processes without a container

| Process | Module | Notes |
|---|---|---|
| `SPLIT_BED` | `modules/local/calling_vardict/split_bed/main.nf` | Only uses coreutils `split` (plus `wc`, `cat`). It declares **no `container` and no `conda` directive**, so it runs in the default execution environment (the Nextflow launcher image/environment). Consider assigning it a lightweight container (e.g. `biocontainers/coreutils`) for full reproducibility. |

---

## Summary of containers with unknown / unconfirmed origin

These are the images we could not fully trace; they should be either documented
or rebuilt from a recipe:

1. `docker.io/bbglab/pysam-0.23.3:latest` — no recipe found; possibly related to the `mutsperpos` recipe.
2. `docker.io/ferriolcalvet/pysam` — untagged personal build; likely the `mutsperpos` recipe.
3. `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5` — personal build; closest recipe is `sigprofiler_matrix`.
4. `docker.io/ferriolcalvet/fgumi:v-retag` — no recipe found anywhere.
5. `biocontainers/seaborn:0.12.2_cv1` — custom-variant BioContainers image; related to the `filtermuts` recipe but build not documented.

