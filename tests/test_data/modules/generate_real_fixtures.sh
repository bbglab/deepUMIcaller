#!/usr/bin/env bash
# Generates small, real (non-synthetic) fixtures for the module-level nf-test
# unit tests under tests/modules/, so those tests can run WITHOUT -stub and
# use `snapshot()` to catch real regressions (e.g. fgbio vs fgumi divergence),
# not just channel/wiring checks.
#
# Run this ON A HOST THAT HAS ACCESS to the source BAMs. By default it points
# at the same reference run used by tests/test_data/input/*.csv (an old run
# parked in fcalvet's scratch space). If you'd rather use a FRESH end-to-end
# run of your own (recommended if that scratch data may be stale/removed, or
# you want fixtures from a specific pipeline version):
#   1. Run the "normal" E2E test with the 5 intermediate BAMs force-published,
#      using tests/test_data/modules/publish_intermediates.config:
#        nextflow run main.nf -profile test,singularity \
#            -c tests/nextflow.config -c tests/test_data/modules/publish_intermediates.config \
#            --input tests/test_data/input/input_test.csv \
#            --outdir /path/to/e2e_normal_out \
#            -work-dir /path/to/e2e_normal_work -resume
#   2. Point this script at it:
#        SRC_BASE=/path/to/e2e_normal_out/bams ./tests/test_data/modules/generate_real_fixtures.sh
#
# It reuses sample B5, cut down to a small subset of real reads (by genomic
# region where the BAM is coordinate-sorted, otherwise by truncating to the
# first few thousand records - see comments below) so the resulting files
# stay small and fast. See docs/test_data.md.
#
# Usage:
#   ./tests/test_data/modules/generate_real_fixtures.sh
#   SRC_BASE=... TARGET_BED=... REF_FASTA=... ./tests/test_data/modules/generate_real_fixtures.sh
#
# Requires: samtools on PATH.

set -euo pipefail

SAMPLE="B5"
SRC_BASE="${SRC_BASE:-/data/bbg/nobackup2/scratch/fcalvet/fetch_duplex_fastqs/bams}"
TARGET_BED="${TARGET_BED:-/data/bbg/nobackup2/scratch/fcalvet/fetch_duplex_fastqs/samplesheet/TP53_target.bed4.bed}"
REF_FASTA="${REF_FASTA:-/data/bbg/datasets/genomes/GRCh38/clean_n_fixed_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.masked.fna}"
OUT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/real"

mkdir -p "$OUT_DIR"

# Resolves to the single *.bam file for $SAMPLE in a directory, without
# assuming an exact name - e.g. SAMTOOLS_SORT (used by all the SORTBAM*
# aliases) renames its output from ".sorted" to ".resorted" whenever the
# input filename already contains ".sorted.", to avoid a name collision. So
# depending on the pipeline path taken, a given stage may be named either way.
find_bam () {
    local dir="$1"
    local matches
    matches=$(find "${dir}" -maxdepth 1 -name "${SAMPLE}*.bam" 2>/dev/null)
    local count
    count=$(echo "${matches}" | grep -c . || true)
    if [ "${count}" -eq 0 ]; then
        echo "ERROR: no ${SAMPLE}*.bam found in ${dir}" >&2
        exit 1
    elif [ "${count}" -gt 1 ]; then
        echo "ERROR: multiple ${SAMPLE}*.bam found in ${dir}, expected exactly one:" >&2
        echo "${matches}" >&2
        exit 1
    fi
    echo "${matches}"
}

# We deliberately do NOT slice by genomic region (samtools view <region>):
# besides only working on coordinate-sorted+indexed BAMs (3 of the 5 stages
# aren't - see conf/modules.config: SORTBAMCLEAN is --template-coordinate,
# CALLDUPLEXCONSENSUSREADS is fgbio's raw unaligned consensus output,
# SORTBAMAMFILTERED is `-n` name-sorted), a single narrow target-panel window
# can genuinely have zero real-world coverage at a given pipeline stage
# (e.g. after duplex consensus collapsing) even when the BAM overall has
# plenty of reads. So instead we just take the header plus the first N
# alignment records everywhere - always non-empty as long as the source BAM
# is, real, and fast (no index needed to read from the start of a file).
N_RECORDS=4000

head_bam () {
    local dir="$1" dst="$2" index_type="${3:-}"
    local src
    src="$(find_bam "${dir}")"
    echo "Truncating (first ${N_RECORDS} records) ${src} -> ${dst}"
    # `head` closing the pipe early would normally SIGPIPE `samtools view`,
    # which combined with `set -o pipefail` aborts the whole script. Drain
    # the rest of the stream with `cat >/dev/null` so samtools always exits
    # cleanly, regardless of how early head stops reading.
    {
        samtools view -H "${src}"
        samtools view "${src}" | { head -n "${N_RECORDS}"; cat >/dev/null; }
    } | samtools view -b - > "${dst}"
    if [ "${index_type}" = "csi" ]; then
        # Only valid because SORTBAMALLMOLECULES/SORTBAMDUPLEXCONS are
        # genome coordinate-sorted, so a prefix of the stream is still
        # validly sorted and can be indexed.
        samtools index -c "${dst}"
    elif [ "${index_type}" = "bai" ]; then
        samtools index "${dst}"
    fi
}

# groupreadsbyumi entry point boundary (FGBIO_GROUPREADSBYUMI input)
head_bam "${SRC_BASE}/SORTBAMCLEAN" \
         "${OUT_DIR}/${SAMPLE}.sortbamclean.bam"

# unmapped_consensus entry point boundary (UNMAP_BAM input)
head_bam "${SRC_BASE}/CALLDUPLEXCONSENSUSREADS" \
         "${OUT_DIR}/${SAMPLE}.consensus.bam"

# filterconsensus entry point boundary (MERGEBAM / FGBIO_FILTERCONSENSUSREADS input)
head_bam "${SRC_BASE}/SORTBAMAMFILTERED" \
         "${OUT_DIR}/${SAMPLE}.sortbamamfiltered.bam"

# allmoleculesfile entry point boundary (ASMINUSXS input, needs .csi)
head_bam "${SRC_BASE}/SORTBAMALLMOLECULES" \
         "${OUT_DIR}/${SAMPLE}.allmolecules.bam" csi

# calling entry point boundary (SAMTOOLS_FILTER real-execution fixture, needs .csi)
head_bam "${SRC_BASE}/SORTBAMDUPLEXCONS" \
         "${OUT_DIR}/${SAMPLE}.duplexcons.bam" csi

# Two (overlapping is fine) subsets for MERGEBAM - it only needs two valid,
# distinct BAM files to merge, not an exact complementary split. MERGEBAM's
# own script always re-indexes its output, which requires coordinate-sorted
# input, so re-sort explicitly here rather than trusting the source's
# existing order (which turned out not to always be reliable in practice).
echo "Deriving two parts of ${SAMPLE}.duplexcons.bam for MERGEBAM..."
samtools view -b -s 42.5 "${OUT_DIR}/${SAMPLE}.duplexcons.bam" | samtools sort -o "${OUT_DIR}/${SAMPLE}.part1.bam" -
samtools view -b -s 43.5 "${OUT_DIR}/${SAMPLE}.duplexcons.bam" | samtools sort -o "${OUT_DIR}/${SAMPLE}.part2.bam" -

# Grouped BAM (fgbio GroupReadsByUmi output) - needed as real input for
# FGBIO_CALLDUPLEXCONSENSUSREADS and FGBIO_COLLECTDUPLEXSEQMETRICS, which are
# also direct fgbio -> fgumi substitution points. Not one of the folders in
# the reference scratch layout, so we derive it ourselves from the already-
# sliced B5.sortbamclean.bam, using the same args as conf/modules.config
# (params.groupreadsbyumi_edits / groupreadsbyumi_min_map_q, both currently
# defaulting to 1 and 10 - update here if those defaults ever change).
if command -v fgbio >/dev/null 2>&1; then
    echo "Deriving ${SAMPLE}.grouped.bam via fgbio GroupReadsByUmi..."
    fgbio -Xmx4g --tmp-dir=. GroupReadsByUmi \
        --edits 1 \
        --min-map-q 10 \
        --strategy Paired \
        --input "${OUT_DIR}/${SAMPLE}.sortbamclean.bam" \
        --output "${OUT_DIR}/${SAMPLE}.grouped.bam" \
        --family-size-histogram "${OUT_DIR}/${SAMPLE}.grouped_umi_histogram.txt"
else
    echo "WARNING: fgbio not found on PATH - skipping ${SAMPLE}.grouped.bam." >&2
    echo "         FGBIO_CALLDUPLEXCONSENSUSREADS and FGBIO_COLLECTDUPLEXSEQMETRICS" >&2
    echo "         real_data tests will fail until it's generated (load/install fgbio" >&2
    echo "         and re-run this script, or run it on a host that has it)." >&2
fi

cp "${TARGET_BED}" "${OUT_DIR}/targets.bed"

# Link the FULL reference, not just one chromosome: the head-truncated BAMs
# above aren't restricted to any single contig (that's the whole point of not
# region-slicing them - see head_bam), so a partial reference can crash fgbio
# the moment it hits a read mapped outside whatever was extracted
# (ArrayIndexOutOfBoundsException in FilterConsensusReads's NM/MD tagging).
# Not committed to the repo either way, so there's no size concern in using
# the whole thing.
echo "Linking full reference genome (not restricted to a single chromosome)..."
ln -sf "${REF_FASTA}" "${OUT_DIR}/reference.fa"
if [ -f "${REF_FASTA}.fai" ]; then
    ln -sf "${REF_FASTA}.fai" "${OUT_DIR}/reference.fa.fai"
else
    samtools faidx "${OUT_DIR}/reference.fa"
fi

echo "Done. Fixtures written to ${OUT_DIR}"
ls -la "${OUT_DIR}"
