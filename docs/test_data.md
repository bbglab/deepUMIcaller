# Test Data Setup

## Origin

The test suite uses publicly available duplex sequencing data from:

> **CRISPR-DS: Combining CRISPR Cas9 enrichment and duplex sequencing for unparalleled
> detection of low frequency mutations**
>
> Schmitt MW, Loeb LA, Salk JJ. *Nat Protoc.* (2016)
> https://pmc.ncbi.nlm.nih.gov/articles/PMC6169890/#s4

Specifically, three TP53 CRISPR-DS samples are used: **B5**, **B7**, and **B13**
(SRA accessions `SRX3224115`, `SRX3224121`, and `SRX3224134` respectively).

---

## Downloading the FASTQ Files

FASTQs can be downloaded using [nf-core/fetchngs](https://nf-co.re/fetchngs).

### 1. Create the accession list

Create a file called `id_list.csv` with the following content:

```csv
SRX3224115
SRX3224121
SRX3224134
```

### 2. Run nf-core/fetchngs

```bash
nextflow run nf-core/fetchngs \
    -profile singularity \
    --input id_list.csv \
    --outdir ./path/to/outdir \
    --download_method sratools
```

**Note**: `nf-core/fetchngs` requires 36 GB of disk space and 6 CPUs.

This will download the FASTQs into `./path/to/outdir/fastq/`, producing files named:

```
SRX3224115_SRR6109268_1.fastq.gz   # B5, read 1
SRX3224115_SRR6109268_2.fastq.gz   # B5, read 2
SRX3224121_SRR6109262_1.fastq.gz   # B7, read 1
SRX3224121_SRR6109262_2.fastq.gz   # B7, read 2
SRX3224134_SRR6109249_1.fastq.gz   # B13, read 1
SRX3224134_SRR6109249_2.fastq.gz   # B13, read 2
```

---

## Preparing Split FASTQs (for multi-file tests)

Tests 3 (`multi-file`) and 5 (`multiAll`) simulate technical replicates by using
split versions of the B7 and B13 FASTQs. Each FASTQ is split into two equal
halves (`part_001` and `part_002`) and placed in a `split_fastqs/` subdirectory.

Repeat the steps below for each of the four FASTQs to split (B7 and B13,
both read 1 and read 2).

**Option A – using `seqkit split2`** (requires [`seqkit`](https://bioinf.shenwei.me/seqkit/)):

```bash
mkdir -p split_fastqs
seqkit split2 -p 2 SRX3224121_SRR6109262_1.fastq.gz --out-dir split_fastqs/
# Repeat for _2.fastq.gz and for SRX3224134_SRR6109249_*.fastq.gz
```

**Option B – using standard Linux tools** (no extra dependencies):

```bash
mkdir -p split_fastqs

# Count the number of reads (4 lines per read in FASTQ)
NLINES=$(zcat SRX3224121_SRR6109262_1.fastq.gz | wc -l)
HALF=$(( NLINES / 2 ))
# Round down to the nearest multiple of 4 to keep reads intact
HALF=$(( HALF - HALF % 4 ))

zcat SRX3224121_SRR6109262_1.fastq.gz | head -n $HALF | gzip > split_fastqs/SRX3224121_SRR6109262_1.part_001.fastq.gz
zcat SRX3224121_SRR6109262_1.fastq.gz | tail -n +$(( HALF + 1 )) | gzip > split_fastqs/SRX3224121_SRR6109262_1.part_002.fastq.gz
# Repeat for _2.fastq.gz and for SRX3224134_SRR6109249_*.fastq.gz
```

The expected output files are:

```
split_fastqs/SRX3224121_SRR6109262_1.part_001.fastq.gz
split_fastqs/SRX3224121_SRR6109262_1.part_002.fastq.gz
split_fastqs/SRX3224121_SRR6109262_2.part_001.fastq.gz
split_fastqs/SRX3224121_SRR6109262_2.part_002.fastq.gz
split_fastqs/SRX3224134_SRR6109249_1.part_001.fastq.gz
split_fastqs/SRX3224134_SRR6109249_1.part_002.fastq.gz
split_fastqs/SRX3224134_SRR6109249_2.part_001.fastq.gz
split_fastqs/SRX3224134_SRR6109249_2.part_002.fastq.gz
```

---

## Preparing BAM Files (for intermediate-step tests)

Tests 6–10 (`groupreadsbyumi`, `filterconsensus`, `calling`, `unmapped_consensus`,
`allmoleculesfile`) start from intermediate pipeline steps and require BAM files
as input. These BAMs are obtained from the output of the `normal` test (Test 1).

After running the `normal` test, locate the relevant BAMs inside the Nextflow
`work/` directory.

```
CALLDUPLEXCONSENSUSREADS   # → input for unmapped_consensus
SORTBAMALLMOLECULES        # → input for allmoleculesfile
SORTBAMAMFILTERED          # → input for filterconsensus
SORTBAMCLEAN               # → input for groupreadsbyumi
SORTBAMDUPLEXCONS          # → input for calling
```

Copy the relevant BAM (and `.csi` index where required) files from those
directories and update the paths in the corresponding input CSVs accordingly.

> ℹ️ More detailed instructions on how to obtain the intermediate BAM files 
> will be added. 

---

## Updating the Input CSV Files

The input CSV files in `tests/test_data/input/` contain **absolute paths** that
must be updated to match the location where you downloaded the data.

After downloading, edit each CSV to replace the path prefix with your local
download directory. For example, if you downloaded to `/my/data/fastq/`:

```bash
# Replace paths in all input CSVs
sed -i 's|/data/bbg/nobackup2/scratch/fcalvet/fetch_duplex_fastqs|/my/data/fastq|g' \
    tests/test_data/input/*.csv
```

---

## Real fixtures for module-level unit tests

The process-level unit tests in `tests/modules/` (see [tests/README.md](../tests/README.md))
each have a `real_data`-tagged variant that runs the real tool (not `-stub`) on
a small, real slice of sample B5 and compares the output against a stored
snapshot. This is what actually catches tool-behaviour regressions (e.g. a
fgbio → fgumi swap changing UMI-grouping output), as opposed to the `-stub`
tests, which only check process wiring.

These real fixtures are **not** committed to the repo (they're real sequencing
data, even if small) - generate them yourself on a host with access to the
same source paths used above:

```bash
./tests/test_data/modules/generate_real_fixtures.sh
```

This truncates the B5 BAMs at each pipeline stage (`SORTBAMCLEAN`,
`CALLDUPLEXCONSENSUSREADS`, `SORTBAMAMFILTERED`, `SORTBAMALLMOLECULES`,
`SORTBAMDUPLEXCONS`) to their first ~4000 alignment records, so the resulting
files stay small while keeping real reads. It intentionally does **not**
slice by genomic region: 3 of the 5 stages aren't coordinate-sorted (so
region queries don't even apply), and even for the 2 that are, a single
narrow target-panel window can have zero real coverage at a given stage
(e.g. after duplex consensus collapsing) - truncating from the start of the
file is simpler and always non-empty. By default it reads from the same old
reference run in fcalvet's scratch space used above; override `SRC_BASE`, `TARGET_BED` and
`REF_FASTA` (as env vars) if your paths differ.

### Using a fresh end-to-end run instead

Three of those five processes (`SORTBAMRAWTEMPCOORDINATE` a.k.a.
`SORTBAMCLEAN`, `CALLCONSENSUSREADS` a.k.a. `CALLDUPLEXCONSENSUSREADS`,
`SORTBAMALLMOLECULES`) have `publishDir enabled: false` in
`conf/modules.config` - they're intermediate and normally only live inside
Nextflow's `work/` directory, not under `--outdir`. If you'd rather generate
the real fixtures from your own fresh run (e.g. to test a specific pipeline
version) instead of that old reference run, force-publish them with
[tests/test_data/modules/publish_intermediates.config](../tests/test_data/modules/publish_intermediates.config):

```bash
nextflow run main.nf -profile test,singularity \
    -c tests/nextflow.config -c tests/test_data/modules/publish_intermediates.config \
    --input tests/test_data/input/input_test.csv \
    --outdir /path/to/e2e_normal_out \
    -work-dir /path/to/e2e_normal_work \
    -resume

SRC_BASE=/path/to/e2e_normal_out/bams ./tests/test_data/modules/generate_real_fixtures.sh
```

