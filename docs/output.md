# bbglab/deepUMIcaller: Output

## Introduction

This document describes the output produced by the pipeline.

## Pipeline overview

- [Directory Structure](#directory-structure)
- [Mutations](#mutations)
- [BAM files](#bam-files)
- [Duplex QC metrics](#duplex-qc-metrics)
- [Basic mutations QC](#basic-mutations-qc)
- [Interal processing files](#internal-files)

## Directory Structure

```{console}
{outdir}
├──callingvardictduplex
├──cohortmutsperpos
├──collectduplexseqmetrics
├──collectduplexseqmetricsontarget
├──computedepth
├──coverageglobal
├──createbedmed
├──discardedcoverageglobal
├──discardedcoveragetargeted
├──familymetrics
├──familymetricsontarget
│   ├── <sample>.med.pdf
│   └── metrics_summary.tsv
├──fastqc
├──multiqc
├──mutations_vcf
│   └── <sample>.med.vcf
├──mutsperpos
├──nsxposition
├──pipeline_info
├──qualimapqcallmolecules
├──qualimapqcduplex
├──qualimapqcraw
├──sigprofplot
├──sigprofplotpur
├──sigprofplotpyr
├──sortbamamfiltered
└──sortbamduplexcons
   └── <sample>.bam
```

## Mutations

Mutations VCF

(optional) Variant annotation

### Key role

- Provides a single-sample VCF with all the mutations detected in duplex reads.
- Format contains information on non-duplex reads in addition to the duplex BAM information. Also the mean position in read and the standard deviation of this value.

### Outputs

- mutations_vcf

## BAM files

### Key role

- BAM files are essential for the downstream analysis to compute the sequencing depth at each position.
- Having the `amfiltered` BAM file allows you to merge sequencing data from the same sample coming from different libraries and doing the calling together.

### Outputs

- sortbamduplexcons
- sortbamamfiltered

## Duplex QC metrics

Metrics to assess the quality of the duplex library prep.

### Key role

- MultiQC provides a summary of the several metrics computed in the pipeline and described here below.
- Family metrics plots are key to make decisions on the optimal sequencing based on the sample diversity.
- FASTQC is important for identifying preexisting problems in the sequencing data.
- QualiMap outputs basic coverage, on target percentage, and other metrics of the different aligned BAM files.
- Discarded coverage files provide a view on the regions of the genome where reads are being lost due to poor mapping.
- Compute depth contains the values of depth per position for each sample. (third column of the .tsv files)
- Coverage global includes a summary of the depth per each region in the global_exons_file provided in the pipeline, in this way there is a bit more useful information on the coverage per region.

### Outputs

- collectduplexseqmetrics
- collectduplexseqmetricsontarget

- discardedcoverageglobal
- discardedcoveragetargeted

- familymetrics
- familymetricsontarget
- fastqc

- multiqc

- qualimapqcallmolecules
- qualimapqcduplex
- qualimapqcraw

- computedepth

- coverageglobal

## Basic mutations QC

### Key role

- Quickly assess the quality of the variants called in the duplex reads. Via simple checks such as the plot of position of the variant in th reads and the mutation frequencies per trinucleotide, you can quickly tell if the variants are occurring based on the expected mutational processes or it is enriched in artifacts. You can also compare the variant calls in purine vs pyrimidine sites and see if there is any obvious differences. No apparent differences a this point does not mean that mutations are all clean.

- cohortmutsperpos includes a representation of the enrichment of mutations at the beginning of the reads for each of the samples in the cohort and for each of the single base substitution types. It also includes a couple of tsv files showing which samples failed the QCs and what are the values of the ratio of mean number of mutations per position found in the first 5 (by default) positions of the read vs the rest of the read.

### Outputs

- cohortmutsperpos
- mutsperpos
- sigprofplot
- sigprofplotpur
- sigprofplotpyr

## Internal files

### Key role

- Each of the following files contains specific information of intermediate processings of the data.

### Outputs

- pipeline_info: compiles basic information of the nextflow run

- callingvardictduplex: contains all the outputs of the variant calling done by VarDict ([click here for more info](https://github.com/AstraZeneca-NGS/VarDictJava))

- createbed: each BED file contains a definition of the regions that have been sequenced in that sample.

- nsxposition: each file has 4 columns: chromosome, position, total_depth (including Ns), and N count.

- readjustregions: this directory contains BED files per sample with information on the regions surrounding mutation calls, targetted regions and others that are used internally for searching and recounting specific mutations. We include it as part of the output just in case it is useful.
