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

## Input and configuration

## Mutations

Mutations VCF
(optional) Variant annotation

### Key role

-

### Outputs

- mutations_vcf

## BAM files

### Key role

-

### Outputs

- sortbamduplexconslow
- sortbamduplexconsmed
- sortbamduplexfiltered

## QC metrics

### Key role

-

### Outputs

- qualimapqc
- qualimapqcduplex
- qualimapqclow
- qualimapqcmed

- discardedcoverageglobal
- discardedcoveragetargeted

- collectduplexseqmetrics
- collectduplexseqmetricsontarget
- familymetrics
- familymetricsontarget
- fastqc

- multiqc

## Basic mutations QC

### Key role

-

### Outputs

- mutsperpos
- sigprofplotlow
- sigprofplotlowpur
- sigprofplotlowpyr
- sigprofplotmedpur
- sigprofplotmedpyr
- sigprofplotmed

Optional:

-

.

## Internal files

### Key role

-

### Outputs

- pipeline_info
- callingvardictlow
- callingvardictmed
- computedepthlow
- computedepthmed
- createbedlow
- createbedmed
- nsxposition
