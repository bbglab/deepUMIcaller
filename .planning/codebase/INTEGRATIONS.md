---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# External Integrations

**Analysis Date:** 2026-10-10

## APIs & External Services

**None at runtime.** The pipeline is fully self-contained: it runs on HPC/local infrastructure with no REST API calls, cloud SDKs, or external service clients in process scripts. All "integrations" are container registries, module repositories, and reference data sources resolved at pipeline *launch* time.

**Container registries (image pulls at task launch):**
- `quay.io` - Default Docker registry for the `docker`/`podman` profiles (`nextflow.config`: `docker.registry = 'quay.io'`); nf-core modules use `quay.io/biocontainers/*` and `quay.io/nf-core/*` images
- `docker.io` - In-house bbglab images (`bbglab/deepcsa-core`, `bbglab/pysam-0.23.3`, `ferriolcalvet/pysam`, `ferriolcalvet/fgumi`, `ferriolcalvet/sigprofilermatrixgenerator`) and `biocontainers/*` images
- `community.wave.seqera.io` - Seqera Wave-built fused container for `modules/local/align_bam/main.nf` (`bwa-mem3_fgumi:ab0c353b848ba934`)
- `depot.galaxyproject.org` - Galaxy Project Singularity mirror used in singularity-mode fallbacks (e.g., `modules/local/samplesheet_check.nf`, `modules/local/mutsperpos/compute/main.nf`)

**Git repositories (module source, resolved by Nextflow at launch):**
- `https://github.com/nf-core/modules.git` - Source of all `modules/nf-core/*` modules; pinned by git SHA in `modules.json` (branch `master`)

**nf-core institutional configs (fetched over HTTPS):**
- `https://raw.githubusercontent.com/nf-core/configs/master` - `params.custom_config_base` in `nextflow.config`; used for `--custom_config_version` institutional config loading

## Data Storage

**Databases:**
- None. No SQL/NoSQL databases. All data is file-based (FASTQ, BAM/CRAM, VCF, BED, TSV) in the Nextflow work directory and published to `params.outdir`.

**Reference data (file-based, local filesystem or shared HPC storage):**
- **iGenomes** - `conf/igenomes.config` builds reference paths from `params.igenomes_base` (BWA/Bowtie2/STAR/Bismark indexes, FASTA, GTF, BED12 for GRCh37, GRCh38, GRCm38, TAIR10, EB2, UMD3.1, WBcel235, CanFam3.1, etc.). In practice the pipeline primarily uses a custom FASTA via `--ref_fasta` (IRB default: `/data/bbg/datasets/genomes/GRCh38/clean_n_fixed_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.masked.fna` in `conf/general_files_IRB.config`).
- **Ensembl VEP cache** - Offline cache directory via `params.vep_cache` (default `null`; IRB default `/data/bbg/datasets/vep` in `conf/general_files_IRB.config`). Defaults: `vep_genome GRCh38`, `vep_species homo_sapiens`, `vep_cache_version 111`, `vep_params "--no_stats --cache --offline --symbol --protein --canonical"` (`nextflow.config`). Consumed by `subworkflows/nf-core/vcf_annotate_ensemblvep/main.nf` → `modules/nf-core/ensemblvep/vep/main.nf`.
- **SigProfiler reference databases** - SigProfilerMatrixGenerator (v1.3.5, container `docker.io/ferriolcalvet/sigprofilermatrixgenerator:1.3.5` in `modules/local/sigprofiler/matrixgenerator/main.nf`) bundles/uses genome-specific signature databases internally; genome assembly passed via `task.ext.genome_assembly` (default `GRCh38`). No external download step in the pipeline.
- **Masking/filter BED files** - Low-complexity, low-mappability, NanoSeq SNP/noise masks, and global exons BEDs supplied as params (IRB paths in `conf/general_files_IRB.config`; blacklist BEDs bundled under `assets/blacklists/` referenced from `conf/igenomes.config`).

**File Storage:**
- Local/POSIX filesystem only (HPC shared storage). Results published under `params.outdir` with routing defined in `conf/modules.config` and `conf/results_outputs.config` (`${params.outdir}/processing_files/<tool>/`, `pipeline_info/`, `metrics/`, `reports/`).
- Singularity image cache: `singularity.cacheDir` / `singularity.libraryDir` = `/data/bbg/datasets/pipelines/nextflow_containers` (`conf/general_files_IRB.config`).

**Caching:**
- Nextflow task-level caching (standard resume); Singularity image cacheDir above. No external cache service.

## Authentication & Identity

**Auth Provider:**
- None. No authentication anywhere in the pipeline. No API keys, tokens, or credentials are read (no `.env` usage; `bin/` scripts take only file/parameter arguments).

## Monitoring & Observability

**Error Tracking:**
- None (no Sentry/external error tracking).

**Logs & Reports (all local files under `${params.outdir}/pipeline_info/`):**
- Nextflow standard execution reports enabled in `nextflow.config`: `timeline`, `report`, `trace`, `dag` (timestamped HTML/TXT files)
- **MultiQC** - Aggregated QC report (`modules/nf-core/multiqc/`), custom config `assets/multiqc_config.yml` (section ordering, `extra_fn_clean_exts` for `.duplex`/`.raw`/`.all_molecules`/`.final` suffixes, `export_plots: true`)
- **Seqera Platform (Tower)** - `tower.yml` at repo root declares which published files are displayed as reports in the Seqera Platform UI (MultiQC report, family metrics plots/TSVs, FastQC, SigProfiler plots, muts-per-position plots). This is a *display* integration; no Tower SDK or `wf-tower` calls exist in the code — launch-time connection is handled by the `tw launch`/Tower agent outside the repo.
- Software versions captured per-process into `versions.yml` (nf-core convention; filtered from publishDir in `conf/modules.config`) and via `custom/dumpsoftwareversions` module (`modules/nf-core/custom/`).

## Notifications

**Email (completion/failure):**
- Implemented in `subworkflows/nf-core/utils_nfcore_pipeline/main.nf` (`completionEmail`): renders the Groovy SimpleTemplate `assets/sendmail_template.txt` and sends via the local `sendmail -t` binary, falling back to `mail` (see lines ~294–316 of that file).
- Wired up (2026-10-21): `main.nf` now invokes `PIPELINE_COMPLETION` after `DEEPUMICALLER`, passing `params.email`, `params.email_on_fail`, `params.plaintext_email`, `params.outdir`, `params.monochrome_logs`, and `DEEPUMICALLER.out.multiqc_report` (new `emit:` on the main workflow). `hook_url` is passed as `null` — webhook notifications are intentionally not wired.
- Triggered when `params.email` or `params.email_on_fail` is set (`subworkflows/local/utils_nfcore_deepumicaller/main.nf` `PIPELINE_COMPLETION` onComplete); HTML by default, plain text with `--plaintext_email`; MultiQC report attached if smaller than `params.max_multiqc_email_size` (default `25.MB`).
- Templates: `assets/email_template.txt` / `assets/email_template.html` (body text/HTML) and `assets/sendmail_template.txt` (MIME wrapper embedding `assets/bbglab_logo.png` as an inline base64 image — replaced the fastquorum-era logo reference on 2026-10-21).
- Requires a functioning local MTA (`sendmail`/`mail`) on the execution host; no SMTP credentials are configured in the repo. `/usr/sbin/sendmail` exists on irbccn26, but whether it relays to external addresses is unverified — test with a real run using `--email`.
- **Best-effort delivery (hardened 2026-10-21):** e-mail problems can never fail the run. `sendMailMessage` (in `subworkflows/nf-core/utils_nfcore_pipeline/main.nf`) checks the send process exit status (a misconfigured MTA no longer logs a false "Sent" message), kills a hung MTA after a 60 s timeout, and surfaces stderr in the warning. Failures fall back sendmail → `mail` (HTML) and are finally caught with a `log.warn` at both the `completionEmail` level and the `PIPELINE_COMPLETION` onComplete handler, so `completionSummary` always runs. Also fixed: `--plaintext_email` was a no-op (exception was constructed but never thrown) and a stray `.sendmail_tmp.html` temp file was written but never used.

**Webhooks/IM (Slack, MS Teams, etc.):**
- `imNotification` is implemented in `subworkflows/nf-core/utils_nfcore_pipeline/main.nf` (Slack/Teams JSON POST), but `main.nf` passes `null` for `hook_url` when invoking `PIPELINE_COMPLETION` — webhook notifications are intentionally not wired (2026-10-21).

## CI/CD & Deployment

**Hosting:**
- Source: GitHub (`https://github.com/bbglab/deepUMIcaller`, per `modules.json` `homePage` and git remote)
- Execution: IRB Barcelona HPC cluster (profile `irbcluster` → `conf/general_files_IRB.config`) or generic local execution (`conf/local.config`)

**CI Pipeline:**
- No GitHub Actions workflows present (`.github/` contains only `ISSUE_TEMPLATE`). Testing is manual/cluster-based via nf-test (`nf-test.config`, `tests/main.nf.test`), with the Singularity profile and a cluster work dir (`DEEPUMICALLER_TEST_WORKDIR`, e.g. `/scratch/bbg/work/deepUMIcaller/nf-tests`).
- No deployment step: the pipeline is consumed by cloning the repo and running `nextflow run bbglab/deepUMIcaller`.

## Environment Configuration

**Required env vars:**
- None strictly required. Optional:
  - `DEEPUMICALLER_TEST_WORKDIR` - overrides nf-test work directory (`nf-test.config`)
  - `PYTHONNOUSERSITE=1`, `R_PROFILE_USER=/.Rprofile`, `R_ENVIRON_USER=/.Renviron` - exported by `nextflow.config` to isolate container Python/R from host libraries

**Key params (all CLI/config, no env):**
- `--input` (samplesheet CSV, validated against `assets/schema_input.json`), `--outdir`, `--step` (`mapping` default; drives pipeline stages), `--ref_fasta`, `--targetsfile`, `--vep_cache` + `--vep_genome/--vep_species/--vep_cache_version`, `--igenomes_base`, `--email/--email_on_fail`

**Secrets location:**
- None in the repository. No `.env` files, no credential files, no tokens. Email delivery relies on the host MTA; container pulls rely on anonymous registry access (or pre-pulled Singularity cache at IRB).

## Webhooks & Callbacks

**Incoming:**
- None.

**Outgoing:**
- Completion email via the local `sendmail` binary when `--email`/`--email_on_fail` is set (see Notifications above). No other outbound webhooks/callbacks; `hook_url` is intentionally unwired.

---

*Integration audit: 2026-10-10*
