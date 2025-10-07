# bbglab/deepUMIcaller: Output

## Introduction

This document collects specific details and comments on tools used in the pipeline or decisions made along the process.

<!-- ## Document overview -->

<!-- - [Basic mutations QC](#basic-mutations-qc)
- [Interal processing files](#internal-files) -->

## VarDict

VarDict is not calling any variants in the first and last 5 bps of a read.

This should be taken into account when setting the clipping parameters, but also note that since we do the pileup for redefining the frequency of each mutation (ALT_DEPTH) the mutations present in these first positions would still be taken into account.
