# deepUMIcaller Input Scenarios Guide

## Overview

The deepUMIcaller pipeline supports multiple input configurations to handle different experimental designs and data structures. This guide explains the various ways to provide input data and configure the pipeline for optimal processing.

## Input File Format

All input configurations use a CSV file with specific columns depending on the pipeline entry point:

### For Standard Mapping (Starting from FASTQ files)

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Sample identifier (used for naming outputs) |
| `fastq_1` | Yes | Path to forward read FASTQ file |
| `fastq_2` | Yes | Path to reverse read FASTQ file |  
| `read_structure` | Yes | UMI read structure specification |
| `parent_dna` | Optional | Patient/individual identifier for multi-sample merging |

### For Intermediate Step Processing

| Entry Point | Required Columns | File Types |  
|-------------|-----------------|------------|  
| `groupreadsbyumi` | `sample`, `bam` | Coordinate-sorted aligned BAM files processed by fgbio GroupByUMI |  
| `unmapped_consensus` | `sample`, `bam` | BAM files with consensus reads that will be realigned |  
| `allmoleculesfile` | `sample`, `duplexbam`, `csi` | BAM with aligned consensus reads missing AS-XS filtering |  
| `filterconsensus` | `sample`, `bam` | BAM with aligned consensus reads only missing a duplex quality filter and the calling |  
| `calling` | `sample`, `duplexbam`, `csi` | Final consensus BAM files + index |

#### Internal outputs that can be used as inputs

You can restart the pipeline from intermediate steps using files produced internally by deepUMIcaller. The table below maps each entry point to the internal process that generates a compatible file, and whether that file is copied to the results folder. This is especially useful when you want to adjust parameters and rerun downstream stages without repeating the most compute‑intensive initial alignments and preprocessing.

| Entry Point | Use deepUMIcaller internal output of | Published to results? |
|-------------|--------------------------------------|------------------------|
| `groupreadsbyumi` | `SORTBAMRAWTEMPCOORDINATE` (coordinate-sorted BAM before UMI grouping) | No (only in work/) |
| `unmapped_consensus` | `CALLCONSENSUSREADS` (consensus BAM prior to alignment/realignment) | No (only in work/) |
| `filterconsensus` | `SORTBAMAMFILTERED` (name-sorted AM-filtered BAM) | Yes → `{outdir}/processing_files/sortbamamfiltered/` |
| `calling` | `SORTBAMDUPLEXCONS` (coordinate-sorted duplex BAM + .csi) | Yes → `{outdir}/duplex_reads_bam/` |
| `allmoleculesfile` | `SORTBAMALLMOLECULES` (coordinate-sorted all-molecules BAM + .csi) | No (only in work/) |

Notes:

- For `calling` and `allmoleculesfile`, provide both the BAM and its CSI index.

- For `filterconsensus`, the pipeline also writes a helper CSV with absolute BAM paths in `{outdir}/pipeline_info/samplesheet_bam_filtered_inputs.csv` that you can pass directly with `--step filterconsensus`.

## Input Scenarios

### 1. Simple scenario

**Use Case**: Standard processing with one FASTQ pair per sample.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure
Sample_A,/path/to/Sample_A_R1.fastq.gz,/path/to/Sample_A_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_R1.fastq.gz,/path/to/Sample_B_R2.fastq.gz,10M1S+T 10M1S+T
Sample_C,/path/to/Sample_C_R1.fastq.gz,/path/to/Sample_C_R2.fastq.gz,10M1S+T 10M1S+T
```

**Pipeline Parameters:**

No additional tuning of parameters required.

**Expected Output:**

- One VCF file per sample: `Sample_A.vcf`, `Sample_B.vcf`, `Sample_C.vcf`
- Independent processing of each sample

---

### 2. Multiple FASTQs for the same sample (same duplex sequencing library)

**Use Case**: When the duplex library preparation of a given sample was sequenced across multiple lanes or runs.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure
Sample_A,/path/to/Sample_A_L1_R1.fastq.gz,/path/to/Sample_A_L1_R2.fastq.gz,10M1S+T 10M1S+T
Sample_A,/path/to/Sample_A_L2_R1.fastq.gz,/path/to/Sample_A_L2_R2.fastq.gz,10M1S+T 10M1S+T
Sample_A,/path/to/Sample_A_L3_R1.fastq.gz,/path/to/Sample_A_L3_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_L1_R1.fastq.gz,/path/to/Sample_B_L1_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_L2_R1.fastq.gz,/path/to/Sample_B_L2_R2.fastq.gz,10M1S+T 10M1S+T
```

**Expected Output:**

- A single VCF per sample: `Sample_A.vcf`, `Sample_B.vcf`
- The files from the different lanes are combined **before** building the duplex reads.
- Independent processing of each sample

---

### 3. Multiple duplex library preparations for the same sample

**Use Case**: When two independent duplex library preparations were done starting from a pool of DNA of the same sample. (e.g., 2x 500 ng library preps to increase duplex coverage).

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure,parent_dna
Sample_A_DNA1_P1,/path/to/Sample_A_DNA1_P1_R1.fastq.gz,/path/to/Sample_A_DNA1_P1_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_A_DNA2_P1,/path/to/Sample_A_DNA2_P1_R1.fastq.gz,/path/to/Sample_A_DNA2_P1_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_B_DNA1_P2,/path/to/Sample_B_DNA1_P2_R1.fastq.gz,/path/to/Sample_B_DNA1_P2_R2.fastq.gz,8M1S+T 8M1S+T,Sample_B_P2
Sample_B_DNA2_P2,/path/to/Sample_B_DNA2_P2_R1.fastq.gz,/path/to/Sample_B_DNA2_P2_R2.fastq.gz,8M1S+T 8M1S+T,Sample_B_P2
```

**Pipeline Parameters:**

If this is the case, note that there is an extra column in the input samplesheet, but apart from this extra column there is no other parameter that should be updated.

No additional tuning of parameters required.

**Expected Output:**

- The files from the different "samples" of the same original DNA, are analyzed independently and are also combined **after** building the duplex reads.
- Individual sample VCFs: `Sample_A_DNA1_P1.vcf`, `Sample_B_DNA1_P2.vcf`, etc.
- Merged sample VCFs: `Sample_A_P1.vcf`, `Sample_B_P2.vcf`

**Benefits:**

- Allows the automated pooling of the sequencing reads in an accurate and robust manner.

---

### 4. Combine scenarios 2 & 3: Multiple lanes, and multiple libraries

**Use Case**: A sample was split across multiple lanes & there are several libraries belonging to the same sample.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,parent_dna
Sample_A_DNA1_P1,/path/to/Sample_A_DNA1_P1_L1_R1.fastq.gz,/path/to/Sample_A_DNA1_P1_L1_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_A_DNA1_P1,/path/to/Sample_A_DNA1_P1_L2_R1.fastq.gz,/path/to/Sample_A_DNA1_P1_L2_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_A_DNA2_P1,/path/to/Sample_A_DNA2_P1_L3_R1.fastq.gz,/path/to/Sample_A_DNA2_P1_L3_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_A_DNA2_P1,/path/to/Sample_A_DNA2_P1_L4_R1.fastq.gz,/path/to/Sample_A_DNA2_P1_L4_R2.fastq.gz,8M1S+T 8M1S+T,Sample_A_P1
Sample_B_DNA1_P2,/path/to/Sample_B_DNA1_P2_R1.fastq.gz,/path/to/Sample_B_DNA1_P2_R2.fastq.gz,8M1S+T 8M1S+T,Sample_B_P2
Sample_B_DNA2_P2,/path/to/Sample_B_DNA2_P2_R1.fastq.gz,/path/to/Sample_B_DNA2_P2_R2.fastq.gz,8M1S+T 8M1S+T,Sample_B_P2
```

**Benefits:**

- Allows the automated pooling of the sequencing reads in an accurate and robust manner.
- Adapts to diverse input scenarios.

**Processing Flow:**

1. **Splitted sequencing lane merging**: before building duplex consensus.
2. **Sample merging**: When appropriate, merge the duplex reads of the multiple libraries of the same original sample.

## Optimization opportunity

Duplex data is known to require to computationally heavy processing.

In certain scenarios, where the targeted panel is very big (i.e. exome) and/or the desired sequencing depth is also very big, the execution of deepUMIcaller might be slow or even get stuck failing to reach the last steps.

Here we list two possibilities for speeding up the processing of the duplex reads that should be used in case of time or memory related issues.

### Chromosome Splitting

**Use Case**: Large datasets requiring computational efficiency.

**Any Input Structure** + Performance Parameter:

```bash
  --split_by_chrom true
```

**Benefits:**

- Parallel processing by chromosome
- Reduced memory requirements per job
- Faster overall execution time
- Same final results as non-split processing

**Compatible with all input scenarios above.**

---

### Split input FASTQs

If the input FASTQs are very big, this will result in long execution times of the first read preprocessing steps. Having the input split across multiple FASTQ files for the same sample will speed up the execution of this steps, or help solve resource-related issues.

We are working to allow the user split the FASTQs within deepUMIcaller, but as of now this is not an option.

## Intermediate Step Entry Points

The deepUMIcaller pipeline supports starting from intermediate processing steps, allowing for:

- **Integration with external pipelines** that produce compatible data formats
- **Leveraging existing processed data** from other UMI or consensus calling tools
- **Component-specific analysis** using deepUMIcaller's specialized modules
- **Workflow modularity** by utilizing specific pipeline components
- **Quality control** and method comparison at different processing stages

### 6. GroupReadsByUmi Step Entry

**Use Case**: Starting from aligned BAM files produced by external alignment pipelines or when you want to apply deepUMIcaller's UMI grouping and consensus calling to existing alignments.

**Input Structure:**

```csv
sample,bam
Sample_A,/path/to/Sample_A.sorted.bam
Sample_B,/path/to/Sample_B.sorted.bam
Sample_C,/path/to/Sample_C.sorted.bam
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input_groupreadsbyumi.csv \
  --outdir results/ \
  --step groupreadsbyumi
```

**Requirements:**

- BAM files must be coordinate-sorted
- BAM files must contain UMI information in read names/tags
- Files must be accessible from compute nodes

**Expected Output:**

- UMI grouping and consensus calling results using deepUMIcaller's algorithms
- Final VCF files with variant calls
- All downstream processing from UMI grouping onward

**Common Scenarios:**

- Using BAM files from external alignment pipelines (BWA, Bowtie2, etc.)
- Comparing UMI grouping methods by applying deepUMIcaller to existing data
- Processing legacy datasets with deepUMIcaller's consensus algorithms

---

### 7. UnmappedConsensus Step Entry

**Use Case**: Starting from BAM files with consensus reads that you want to realign to you genome of interest. The input BAM file can either be unmapped or mapped to another genome assembly.

**Input Structure:**

```csv
sample,bam
Sample_A,/path/to/Sample_A.unmapped.bam
Sample_B,/path/to/Sample_B.unmapped.bam
Sample_C,/path/to/Sample_C.unmapped.bam
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input_unmapped_consensus.csv \
  --outdir results/ \
  --step unmapped_consensus
```

**Requirements:**

- BAM files must contain unmapped reads with UMI information
- Files must be properly formatted for consensus calling
- UMI tags must be present in read names or BAM tags

**Expected Output:**

- Unmapped consensus reads processed using deepUMIcaller's algorithms
- Downstream processing from consensus generation onward
- Final variant calling results

**Common Scenarios:**

- Processing unmapped reads from external pipelines
- Handling reads that failed initial alignment but contain valid UMI information
- Quality control on unmapped read consensus generation

---

### 8. AllMoleculesFile Step Entry

**Use Case**: Starting from final processed duplex BAM files to generate comprehensive all-molecules output files, useful for detailed molecular analysis and quality metrics across all UMI families.

**Input Structure:**

```csv
sample,duplexbam,csi
Sample_A,/path/to/Sample_A.duplex.bam,/path/to/Sample_A.duplex.bam.csi
Sample_B,/path/to/Sample_B.duplex.bam,/path/to/Sample_B.duplex.bam.csi
Sample_C,/path/to/Sample_C.duplex.bam,/path/to/Sample_C.duplex.bam.csi
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input_allmoleculesfile.csv \
  --outdir results/ \
  --step allmoleculesfile
```

**Requirements:**

- BAM files must contain fully processed duplex consensus reads
- CSI index files must be provided for each BAM
- Files must be ready for all-molecules analysis

**Expected Output:**

- Comprehensive all-molecules output files
- Detailed metrics per UMI family
- Quality control statistics across all molecular families
- Analysis-ready data for downstream investigations

**Common Scenarios:**

- Generating comprehensive molecular statistics from existing duplex data
- Quality control and validation of UMI family processing
- Detailed analysis of molecular coverage and family sizes
- Research applications requiring per-molecule information

---

### 9. FilterConsensus Step Entry

**Use Case**: Starting from consensus BAM files produced by external UMI processing tools, or when you want to apply deepUMIcaller's consensus filtering and variant calling to existing consensus reads.

**Input Structure:**

```csv
sample,bam
Sample_A,/path/to/Sample_A.consensus.bam
Sample_B,/path/to/Sample_B.consensus.bam
Sample_C,/path/to/Sample_C.consensus.bam
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input_filterconsensus.csv \
  --outdir results/ \
  --step filterconsensus
```

**Requirements:**

- BAM files must contain consensus reads from UMI processing
- Files should be from completed UMI grouping and consensus calling
- Proper consensus read formatting required

**Expected Output:**

- Consensus read filtering and quality control using deepUMIcaller's methods
- Final variant calling results
- High-quality filtered VCF outputs

Tip: deepUMIcaller writes a helper samplesheet at `{outdir}/sortbamamfiltered/samplesheet_bam_filtered_inputs.csv` with absolute BAM paths. You can use this directly with `--step filterconsensus` to restart from these AM-filtered BAMs.

**Common Scenarios:**

- Processing consensus BAMs from external UMI tools (UMI-tools, CGAT-core, etc.)
- Applying deepUMIcaller's filtering algorithms to existing consensus data
- Benchmarking variant calling performance on standardized consensus inputs

---

### 10. Calling Step Entry

**Use Case**: Starting from final processed duplex BAM files produced by external consensus pipelines, or when you specifically want to use deepUMIcaller's variant calling algorithms on high-quality duplex consensus data.

**Input Structure:**

```csv
sample,duplexbam,csi
Sample_A,/path/to/Sample_A.duplex.bam,/path/to/Sample_A.duplex.bam.csi
Sample_B,/path/to/Sample_B.duplex.bam,/path/to/Sample_B.duplex.bam.csi
Sample_C,/path/to/Sample_C.duplex.bam,/path/to/Sample_C.duplex.bam.csi
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input_calling.csv \
  --outdir results/ \
  --step calling
```

**Requirements:**

- BAM files must contain fully processed duplex consensus reads
- CSI index files must be provided for each BAM
- Files must be ready for variant calling

**Expected Output:**

- Variant calling results using deepUMIcaller's algorithms
- VCF files with detected variants
- Quality metrics and statistics

**Common Scenarios:**

- Using duplex BAMs from external consensus calling pipelines
- Applying deepUMIcaller's variant calling to standardized duplex data
- Comparing variant calling methods using the same input consensus data
- Processing duplex consensus data from specialized UMI workflows

---


## Parameter Interactions

### Key Parameters

| Parameter | Default | Effect |
|-----------|---------|--------|
| `step` | `"mapping"` | Pipeline entry point (mapping, groupreadsbyumi, unmapped_consensus, filterconsensus, calling, allmoleculesfile) |
| `split_by_chrom` | `false` | Enables chromosome-based parallelization |
| `parent_dna` column | - | Enables biological replicate grouping |

### Parameter Combinations

| Step | FASTQs splitted at origin | Multiple libraries per sample | Chromosome Splitting | Use Case |
|------|---------------------------|-------------------------------|----------------------|----------|
| mapping | ❌ | ❌ | ❌ | Basic single-sample processing |
| mapping | ✅ | ❌ | ❌ | 2 |
| mapping | ❌ | ✅ | ❌ | 3 |
| mapping | ✅ | ✅ | ❌ | 4 |
| mapping | ✅/❌ | ✅/❌ | ✅ | Parallelized processing of information per chromosome |
| groupreadsbyumi | ❌ | ❌ | ❌ | UMI grouping restart |
| unmapped_consensus | ❌ | ❌ | ❌ | Unmapped consensus processing |
| filterconsensus | ❌ | ❌ | ❌ | Consensus filtering restart |
| calling | ❌ | ❌ | ❌ | Variant calling only |
| allmoleculesfile | ❌ | ❌ | ❌ | All-molecules analysis |

## Best Practices

### 1. File Path Management

- Use **absolute paths** for FASTQ files
- Ensure all files are accessible from compute nodes
- Consider using symbolic links for large datasets

### 2. Sample Naming

- Use **consistent naming conventions**
- Avoid special characters (dots, spaces, brackets, etc.)
- Keep names descriptive but concise

### 3. Parent DNA Identifiers

- Use **consistent patient/individual IDs**
- Ensure proper grouping of related samples
- Consider using study-specific prefixes

### 4. Performance Considerations

- Enable `split_by_chrom` for heavy samples (>100 GBs of raw reads )
- Use technical replicate merging when appropriate
- Monitor resource usage and adjust accordingly

## Validation and Quality Control

### Input Validation

The pipeline automatically validates:

- CSV file format and required columns
- FASTQ file existence and accessibility
- Sample name uniqueness within processing groups
- Parent DNA identifier consistency

### Quality Metrics

- **Coverage depth** per sample and merged group
- **Variant calling statistics** per processing level
- **Precision metrics** comparing processing strategies

## Troubleshooting

### Common Issues

1. **Duplicate sample names**: Ensure unique sample identifiers
2. **Memory issues**: Consider chromosome splitting for large datasets
3. **Invalid step parameter**: Use valid step names (mapping, groupreadsbyumi, unmapped_consensus, filterconsensus, calling, allmoleculesfile)
4. **Incorrect file format**: Ensure BAM files match the expected processing stage

For additional support or questions about input configurations, consult the main deepUMIcaller documentation or contact the development team.
