---
last_mapped_commit: d3718542cbe52b31362e4dc42c1c753054710b9d
last_mapped_at: 2026-10-10
---
# Testing Patterns

**Analysis Date:** 2026-10-10

## Test Framework

**Runner:**
- **nf-test** (pipeline-level integration testing) — the only test framework in the repo
- Config: `nf-test.config` (repo root)
- No pytest/unittest for Python helpers; no unit test framework for Groovy functions

**nf-test configuration** (`nf-test.config`):

```groovy
config {
    testsDir "tests"
    workDir System.getenv("DEEPUMICALLER_TEST_WORKDIR") ?: ".nf-test"
    configFile "tests/nextflow.config"
    profile "test,singularity"          // activate test + singularity profiles while testing
    ignore 'modules/nf-core/**/*', 'subworkflows/nf-core/**/*'   // nf-core modules have their own CI
    plugins {
        load "nft-utils@0.0.3"
    }
}
```

**Run Commands:**

```bash
nf-test test                          # Run all tests
nf-test test --tag "normal"           # Run one test by tag
nf-test test --tag "multiAll"         # Comprehensive scenario
nf-test test --verbose                # Verbose output
```

Tags available: `normal`, `split_fastq`, `split_by_chrom`, `multi-file`, `multi-sample`, `multiAll`, `groupreadsbyumi`, `filterconsensus`, `calling`, `unmapped_consensus`, `allmoleculesfile` (see `tests/README.md` for the full catalog).

**Workdir override:** set `DEEPUMICALLER_TEST_WORKDIR` to point nf-test work at cluster scratch (recommended on the IRB cluster); defaults to `.nf-test/` in the repo root (gitignored).

## Test File Organization

**Location:**
- Single pipeline-level suite: `tests/main.nf.test` (~1100 lines, 11 tests)
- No per-module nf-test files; `modules/nf-core/**` and `subworkflows/nf-core/**` are excluded from testing via `ignore` in `nf-test.config`

**Naming:**
- Suite file: `tests/main.nf.test`
- Test names follow `TEST <n>. <Category> - <description>` (e.g., `"TEST 4. Technical replicates - multi-file sample processing"`)
- Tags are lowercase kebab/snake: `multi-file`, `split_by_chrom`

**Structure:**

```
tests/
├── main.nf.test           # nf-test pipeline suite (11 tests)
├── nextflow.config        # test-specific Nextflow config (loaded by nf-test.config)
├── README.md              # test suite documentation + tag catalog
├── 2026-10-07.results.csv # historical nf-test run results (PASSED/FAILED per test)
├── 2026-10-09.results.csv # historical nf-test run results
└── test_data/
    ├── input/             # 10 samplesheet CSVs (one per test scenario)
    └── expected_output/   # reference VCFs (B5.vcf, B7.vcf, B13.vcf, Bsample.vcf)
```

## Test Structure

**Suite organization** (`tests/main.nf.test`):

```groovy
nextflow_pipeline {
    name "Test Workflow main.nf"
    script "main.nf"
    config "./nextflow.config"

    test("TEST 1. Basic functionality - standard processing") {
        tag "normal"

        when {
            params {
                input  = "${projectDir}/tests/test_data/input/input_test.csv"
                outdir = "$outputDir"
            }
        }
        then {
            assert workflow.success : "Pipeline should complete without errors"
            // ... VCF precision validation ...
        }
    }
}
```

**Patterns:**
- **Setup:** all inputs are passed via `when { params { ... } }` — no channel setup in tests; each test gets its own `outdir` subdirectory (`${outputDir}/tests_results_<scenario>/`)
- **Assertion:** `assert workflow.success : "message"` first, then file-level assertions with descriptive failure messages prefixed by category: `VALIDATION ERROR:`, `REFERENCE ERROR:`, `SCRIPT ERROR:`, `COMPARISON SCRIPT ERROR:`, `PRECISION FAILURE:`
- **Teardown:** none (nf-test workdirs are disposable)

**Test scenarios covered** (one samplesheet per scenario in `tests/test_data/input/`):
| Test | Tag | Scenario | Input CSV |
|------|-----|----------|-----------|
| 1 | `normal` | Single FASTQ pair per sample | `input_test.csv` |
| 2 | `split_fastq` | Internal FASTQ splitting (`run_splitfastq=true`, `splitfastq_parts=2`) | `input_test.csv` |
| 3 | `split_by_chrom` | Chromosome parallelization (`split_by_chrom=true`) | `input_test_new1.csv` |
| 4 | `multi-file` | Technical replicates (lanes) | `input_test_multi-file.csv` |
| 5 | `multi-sample` | Biological replicates (`parent_dna` grouping) | `input_test_multi-sample.csv` |
| 6 | `multiAll` | All combined | `input_test_multiAll.csv` |
| 7 | `groupreadsbyumi` | Entry at `step=groupreadsbyumi` | `input_groupreadsbyumi.csv` |
| 8 | `filterconsensus` | Entry at `step=filterconsensus` | `input_filterconsensus.csv` |
| 9 | `calling` | Entry at `step=calling` | `input_calling.csv` |
| 10 | `unmapped_consensus` | Entry at `step=unmapped_consensus` | `input_unmapped_consensus.csv` |
| 11 | `allmoleculesfile` | Entry at `step=allmoleculesfile` | `input_allmoleculesfile.csv` |

## VCF Precision Validation (core assertion pattern)

Every test validates output VCFs against reference VCFs by calling a Python comparison script from Groovy and asserting a ≥99% precision threshold:

```groovy
def currentDir = file("${params.outdir}/processing_files/callingvardictduplex")
def expectedDir = file("tests/test_data/expected_output")
def currentVcfs = currentDir.listFiles().findAll { it.name.endsWith('.vcf') }
assert currentVcfs.size() > 0 : "VALIDATION ERROR: No VCF files generated in output directory: ${currentDir}"

currentVcfs.each { currentVcf ->
    def expectedVcf = new File(expectedDir, currentVcf.name)
    assert expectedVcf.exists() : "REFERENCE ERROR: Expected reference file missing: ${expectedVcf.name}"

    def scriptPath = new File("${projectDir}/bin/test_utilities/compare_vcfs_detailed.py")
    def cmd = "python3 ${scriptPath} ${currentVcf} ${expectedVcf}"
    def proc = cmd.execute()
    proc.waitFor()
    def output = proc.in.text.trim()
    assert proc.exitValue() == 0 : "COMPARISON SCRIPT ERROR: ... exit code ${proc.exitValue()}. Error: ${proc.err.text}"
    assert output.isNumber() : "COMPARISON SCRIPT ERROR: Non-numeric output ...: '${output}'"
    def precision = output as Double

    println "\n[TEST 1 - BASIC] VCF precision for ${currentVcf.name}: ${precision}%"
    assert precision >= 99.0 : "PRECISION FAILURE: VCF ${currentVcf.name} precision ${precision}% below 99% threshold"
}
```

(`tests/main.nf.test:73-105`, repeated in all 11 tests)

**Comparison scripts:**
- `bin/test_utilities/compare_vcfs_detailed.py` — computes Jaccard/precision/recall/F1 between two VCFs; default output is a single numeric percentage (parsed by tests); `debug` mode prints variant lists. Variants keyed as `(chrom, pos, ref, alt)` tuples; malformed lines silently skipped.
- `bin/test_utilities/validate_vcf_biology.py` — sanity checks (valid chromosome names, numeric positions, valid REF/ALT bases)
- `bin/compare_vcfs.py` — older minimal Jaccard-only comparator (superseded by `compare_vcfs_detailed.py`)

## Test Data Organization

**Input samplesheets** (`tests/test_data/input/*.csv`) reference real SRA FASTQs on IRB cluster storage:

```csv
sample,fastq_1,fastq_2,read_structure
B5,/data/bbg/nobackup2/scratch/fcalvet/fetch_duplex_fastqs/fastq/SRX3224115_SRR6109268_1.fastq.gz,/data/bbg/nobackup2/scratch/fcalvet/fetch_duplex_fastqs/fastq/SRX3224115_SRR6109268_2.fastq.gz,10M1S+T 10M1S+T
```

(`tests/test_data/input/input_test.csv`)

**Expected outputs** (`tests/test_data/expected_output/`): reference VCFs named by sample (`B5.vcf`, `B7.vcf`, `B13.vcf`, `Bsample.vcf`); matched by filename against produced VCFs. Also contains `fgbio/` subdirectory and `vcf2maf.py` helper.

**Test-specific Nextflow config** (`tests/nextflow.config`):
- Loads cluster references: `includeConfig '../conf/general_files_IRB.config'`
- Overrides test params: `targetsfile`, `right_clip=30`, `asminusxs_thres=10`, `final_read_size=270`, `filter_min_reads_duplex="2,1,1"`
- Executor: `slurm` with `errorStrategy = 'retry'`, `maxRetries = 2`
- Disables publishing for `PATCH_DEPTH`, sets `MUTSPERPOS` to `errorStrategy = 'ignore'`
- Enables trace/timeline/report/DAG under `${params.outdir}/pipeline_info/`
- Singularity cacheDir/libraryDir pointed at `/data/bbg/datasets/pipelines/nextflow_containers`

**Important:** test data paths are cluster-specific (IRB `/data/bbg/...`). Tests cannot run off-cluster without replacing these paths.

## Mocking

**Framework:** none.

**What is mocked:** nothing. This is a bioinformatics pipeline integration suite — real tools run in Singularity containers against real (small) data.

**What NOT to mock:** pipeline execution, tool outputs. Determinism is achieved via fixed reference VCFs and a numeric precision threshold rather than exact snapshot matching.

## Fixtures and Factories

**Test Data:**
- Static CSV samplesheets per scenario (no factories/generators)
- `assets/samplesheet.csv` provides a minimal example for users; `assets/schema_input.json` defines valid columns
- `testing/` directory exists but is empty (gitignored via `.gitignore` `testing*` entries) — used for ad-hoc local runs

**Location:** `tests/test_data/input/` (inputs), `tests/test_data/expected_output/` (references)

## Coverage

**Requirements:** None enforced (no coverage tooling for Nextflow or Python).

**Proxy for coverage:** the 11-test matrix exercises all `step` entry points (`mapping`, `groupreadsbyumi`, `unmapped_consensus`, `filterconsensus`, `calling`, `allmoleculesfile`), all merging modes (multi-file, multi-sample), and both parallelization flags (`split_by_chrom`, `run_splitfastq`).

**View Results:**

```bash
nf-test test --verbose
```

Historical run results are committed as CSV: `tests/2026-10-07.results.csv`, `tests/2026-10-09.results.csv` (columns: filename, testsuite, type, test, result, start, end, time). These show per-test PASSED/FAILED status and wall time (tests take ~20-40 min each on the cluster).

## Test Types

**Pipeline Integration Tests:**
- All testing is end-to-end: `nextflow_pipeline` blocks run `main.nf` with the `test,singularity` profiles
- Runtime is significant (each test 20-40 min; full suite several hours per `tests/2026-10-09.results.csv`)

**Unit Tests:**
- Not used. Python helpers in `bin/` have no unit tests; Groovy helpers in `modules/local/utils.nf` and `subworkflows/local/input_check.nf` (`create_input_channel`) have no unit tests.

**E2E Tests:**
- The nf-test suite *is* the E2E layer.

**Module-level tests:**
- None for `modules/local/`; nf-core modules under `modules/nf-core/` are explicitly excluded (`nf-test.config` `ignore`) because they carry their own upstream CI.

## Common Patterns

**Parameterized scenario testing:**

```groovy
test("TEST 3. Performance optimization - chromosome-based parallelization") {
    tag "split_by_chrom"
    when {
        params {
            input          = "${projectDir}/tests/test_data/input/input_test_new1.csv"
            outdir         = "${outputDir}/tests_results_splitChrom/"
            split_by_chrom = true
        }
    }
    then { /* same precision assertions as TEST 1 */ }
}
```

**Numeric-threshold regression testing** (instead of exact file diffing):

```groovy
assert precision >= 99.0 : "PRECISION FAILURE: ..."
```

This tolerates non-deterministic ordering/metadata differences in VCFs while catching biological regressions.

**Intermediate-step entry testing:** each `step` value gets its own samplesheet + test, validating that restarting mid-pipeline produces equivalent final VCFs.

## Gaps and Notes

- **No CI:** `.github/workflows/` does not exist — tests are run manually on the IRB cluster; results CSVs are committed as records
- **No Python unit tests:** `bin/*.py` (3,476 lines total) is only exercised indirectly through pipeline runs
- **No module-level nf-test files:** local processes are only tested via full-pipeline runs
- **Cluster-bound test data:** inputs/expected outputs reference `/data/bbg/...` paths; portability requires re-hosting data
- **Duplicated assertion blocks:** the ~30-line precision-validation block is copy-pasted across all 11 tests in `tests/main.nf.test` — extract into a shared Groovy closure if adding tests
- **`testing/` directory is empty** and gitignored; do not place permanent test assets there

---

*Testing analysis: 2026-10-10*
