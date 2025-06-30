# bbglab/deepUMIcaller: Output

## Introduction

This document describes the output produced by the pipeline.

## Pipeline overview

- [Directory Structure](#directory-structure)
- [Input and configuration](#input-and-configuration)
- [Mutations](#mutations)
- [BAM files](#bam-files)
- [QC metrics](#qc-metrics)
- [Basic mutations QC](#basic-mutations-qc)
- [Interal processing files](#internal-files)

## Directory Structure

```{console}
{outdir}
├──callingvardictmed
├──collectduplexseqmetrics
├──collectduplexseqmetricsontarget
├──computedepthmed
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
├──qualimapqc
├──qualimapqcduplex
├──qualimapqcmed
├──sigprofplotmed
├──sigprofplotmedpur
├──sigprofplotmedpyr
├──sortbamduplexconsmed
│   └── <sample>.bam
└──sortbamduplexfiltered
```

In case the pipeline was run, activating the `duplex_high_conf` and `duplex_low_conf` parameters, several folders would be replicated with the respective confidence level's data.

## Input and configuration

This information is available in the usage.md document.

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
- Having the `duplexfiltered` BAM file allows you to merge sequencing data from the same sample coming from different libraries and doing the calling together.

### Outputs

- sortbamduplexconsmed
- sortbamduplexfiltered

## QC metrics

Metrics to assess the quality of the duplex library prep.

### Key role

- MultiQC provides a summary of the several metrics computed in the pipeline and described here below.
- Family metrics plots are key to make decisions on the optimal sequencing based on the sample diversity.
- FASTQC is important for identifying preexisting problems in the sequencing data.
- QualiMap outputs basic coverage, on target percentage, and other metrics of the different aligned BAM files.
- Discarded coverage files provide a view on the regions of the genome where reads are being lost due to poor mapping.

### Outputs

- collectduplexseqmetrics
- collectduplexseqmetricsontarget

- discardedcoverageglobal
- discardedcoveragetargeted

- familymetrics
- familymetricsontarget
- fastqc

- multiqc

- qualimapqc
- qualimapqcduplex
- qualimapqcmed

## Basic mutations QC

### Key role

- Quickly assess the quality of the variants called in the duplex reads. Via simple checks such as the plot of position of the variant in th reads and the mutation frequencies per trinucleotide, you can quickly tell if the variants are occurring based on the expected mutational processes or it is enriched in artifacts. You can also compare the variant calls in purine vs pyrimidine sites and see if there is any obvious differences. No apparent differences a this point does not mean that mutations are all clean.

### Outputs

- mutsperpos
- sigprofplotmed
- sigprofplotmedpur
- sigprofplotmedpyr

## Internal files

### Key role

- Each of the following files contains specific information of intermediate processings of the data.

### Outputs

- pipeline_info
- callingvardictmed
- computedepthmed
- createbedmed
- nsxposition
