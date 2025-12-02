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
| `groupreadsbyumi` | `SORTBAMCLEAN` (coordinate-sorted BAM before UMI grouping) | No (only in work/) |
| `unmapped_consensus` | `CALLDUPLEXCONSENSUSREADS` (consensus BAM prior to alignment/realignment) | No (only in work/) |
| `filterconsensus` | `SORTBAMAMFILTERED` (name-sorted AM-filtered BAM) | Yes → `{outdir}/sortbamamfiltered/` |
| `calling` | `SORTBAMDUPLEXCONS` (coordinate-sorted duplex BAM + .csi) | Yes → `{outdir}/sortbamduplexcons/` |
| `allmoleculesfile` | `SORTBAMALLMOLECULES` (coordinate-sorted all-molecules BAM + .csi) | No (only in work/) |

Notes:
- For `calling` and `allmoleculesfile`, provide both the BAM and its CSI index.
- For `filterconsensus`, the pipeline also writes a helper CSV with absolute BAM paths in `{outdir}/sortbamamfiltered/samplesheet_bam_filtered_inputs.csv` that you can pass directly with `--step filterconsensus`.

## Input Scenarios

### 1. Single Sample Processing (Basic)

**Use Case**: Standard processing with one FASTQ pair per sample.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure
Sample_A,/path/to/Sample_A_R1.fastq.gz,/path/to/Sample_A_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_R1.fastq.gz,/path/to/Sample_B_R2.fastq.gz,10M1S+T 10M1S+T
Sample_C,/path/to/Sample_C_R1.fastq.gz,/path/to/Sample_C_R2.fastq.gz,10M1S+T 10M1S+T
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input.csv \
  --outdir results/
```

**Expected Output:**

- One VCF file per sample: `Sample_A.vcf`, `Sample_B.vcf`, `Sample_C.vcf`
- Independent processing of each sample
- No cross-sample merging

---

### 2. Technical Replicates (Multi-lane Merging)

**Use Case**: When the same biological sample was sequenced across multiple lanes or runs.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure
Sample_A,/path/to/Sample_A_L1_R1.fastq.gz,/path/to/Sample_A_L1_R2.fastq.gz,10M1S+T 10M1S+T
Sample_A,/path/to/Sample_A_L2_R1.fastq.gz,/path/to/Sample_A_L2_R2.fastq.gz,10M1S+T 10M1S+T
Sample_A,/path/to/Sample_A_L3_R1.fastq.gz,/path/to/Sample_A_L3_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_L1_R1.fastq.gz,/path/to/Sample_B_L1_R2.fastq.gz,10M1S+T 10M1S+T
Sample_B,/path/to/Sample_B_L2_R1.fastq.gz,/path/to/Sample_B_L2_R2.fastq.gz,10M1S+T 10M1S+T
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input.csv \
  --outdir results/ \
  --splitted_original_sample true
```

**Expected Output:**

- Merged VCF per biological sample: `Sample_A.vcf`, `Sample_B.vcf`
- Improved sensitivity from combined coverage
- Lane-level BAM files merged before variant calling

**Benefits:**

- Increased sequencing depth per sample
- Better variant detection sensitivity
- Reduced false negative rates

---

### 3. Biological Replicates (Multi-sample Patient Processing)

**Use Case**: Multiple samples from the same patient/individual (e.g., different tissues, timepoints).

**Input Structure:**

```csv
sample,fastq_1,fastq_2,read_structure,parent_dna
Sample_A_P1,/path/to/Sample_A_P1_R1.fastq.gz,/path/to/Sample_A_P1_R2.fastq.gz,10M1S+T 10M1S+T,Patient_001
Sample_B_P1,/path/to/Sample_B_P1_R1.fastq.gz,/path/to/Sample_B_P1_R2.fastq.gz,10M1S+T 10M1S+T,Patient_001
Sample_C_P1,/path/to/Sample_C_P1_R1.fastq.gz,/path/to/Sample_C_P1_R2.fastq.gz,10M1S+T 10M1S+T,Patient_001
Sample_A_P2,/path/to/Sample_A_P2_R1.fastq.gz,/path/to/Sample_A_P2_R2.fastq.gz,10M1S+T 10M1S+T,Patient_002
Sample_B_P2,/path/to/Sample_B_P2_R1.fastq.gz,/path/to/Sample_B_P2_R2.fastq.gz,10M1S+T 10M1S+T,Patient_002
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input.csv \
  --outdir results/
```

**Expected Output:**

- Sample-level VCFs: `Sample_A_P1.vcf`, `Sample_B_P1.vcf`, etc.
- Patient-level merged VCFs: `Patient_001.vcf`, `Patient_002.vcf`
- Enhanced variant calling through multi-sample evidence

**Benefits:**

- Patient-level variant consolidation
- Cross-sample variant validation
- Comprehensive genomic profiling per patient

---

### 4. Combined Technical + Biological Replicates

**Use Case**: Complex experimental designs with both technical and biological replicates.

**Input Structure:**

```csv
sample,fastq_1,fastq_2,parent_dna
Sample_A_P1,/path/to/Sample_A_P1_L1_R1.fastq.gz,/path/to/Sample_A_P1_L1_R2.fastq.gz,Patient_001
Sample_A_P1,/path/to/Sample_A_P1_L2_R1.fastq.gz,/path/to/Sample_A_P1_L2_R2.fastq.gz,Patient_001
Sample_B_P1,/path/to/Sample_B_P1_L1_R1.fastq.gz,/path/to/Sample_B_P1_L1_R2.fastq.gz,Patient_001
Sample_B_P1,/path/to/Sample_B_P1_L2_R1.fastq.gz,/path/to/Sample_B_P1_L2_R2.fastq.gz,Patient_001
Sample_A_P2,/path/to/Sample_A_P2_L1_R1.fastq.gz,/path/to/Sample_A_P2_L1_R2.fastq.gz,Patient_002
Sample_A_P2,/path/to/Sample_A_P2_L2_R1.fastq.gz,/path/to/Sample_A_P2_L2_R2.fastq.gz,Patient_002
```

**Pipeline Parameters:**

```bash
nextflow run main.nf \
  --input input.csv \
  --outdir results/ \
  --splitted_original_sample true
```

**Expected Output:**

- Merged sample-level VCFs: `Sample_A_P1.vcf`, `Sample_B_P1.vcf`, etc.
- Patient-level consolidated VCFs: `Patient_001.vcf`, `Patient_002.vcf`
- Maximum sensitivity from both technical and biological aggregation

**Processing Flow:**

1. **Lane Merging**: Technical replicates merged per sample
2. **Sample Processing**: Individual variant calling per merged sample
3. **Patient Consolidation**: Biological samples merged per patient

---

### 5. Performance Optimization (Chromosome Splitting)

**Use Case**: Large datasets requiring computational efficiency.

**Any Input Structure** + Performance Parameter:

```bash
nextflow run main.nf \
  --input input.csv \
  --outdir results/ \
  --split_by_chrom true
```

**Benefits:**

- Parallel processing by chromosome
- Reduced memory requirements per job
- Faster overall execution time
- Same final results as non-split processing

**Compatible with all input scenarios above.**

---

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
| `splitted_original_sample` | `false` | Enables technical replicate merging |
| `split_by_chrom` | `false` | Enables chromosome-based parallelization |
| `parent_dna` column | - | Enables biological replicate grouping |

### Parameter Combinations

| Step | Technical Replicates | Biological Replicates | Chromosome Splitting | Use Case |
|------|---------------------|---------------------|---------------------|----------|
| mapping | ❌ | ❌ | ❌ | Basic single-sample processing |
| mapping | ❌ | ❌ | ✅ | Basic + performance optimization |
| mapping | ✅ | ❌ | ❌ | Technical replicate merging |
| mapping | ❌ | ✅ | ❌ | Biological replicate grouping |
| mapping | ✅ | ✅ | ❌ | Multi-level processing |
| mapping | ✅ | ✅ | ✅ | Comprehensive + optimized |
| groupreadsbyumi | ❌ | ❌ | ❌ | UMI grouping restart |
| unmapped_consensus | ❌ | ❌ | ❌ | Unmapped consensus processing |
| filterconsensus | ❌ | ❌ | ✅ | Consensus filtering restart |
| calling | ❌ | ❌ | ❌ | Variant calling only |
| allmoleculesfile | ❌ | ❌ | ❌ | All-molecules analysis |

## Best Practices

### 1. File Path Management

- Use **absolute paths** for FASTQ files
- Ensure all files are accessible from compute nodes
- Consider using symbolic links for large datasets

### 2. Sample Naming

- Use **consistent naming conventions**
- Avoid special characters (spaces, brackets, etc.)
- Keep names descriptive but concise

### 3. Parent DNA Identifiers

- Use **consistent patient/individual IDs**
- Ensure proper grouping of related samples
- Consider using study-specific prefixes

### 4. Performance Considerations

- Enable `split_by_chrom` for large datasets (>50 samples)
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

### Common Issues:

1. **Duplicate sample names**: Ensure unique sample identifiers
2. **Missing parent_dna**: Required for biological replicate merging
3. **File path errors**: Verify FASTQ/BAM file accessibility
4. **Memory issues**: Consider chromosome splitting for large datasets
5. **Invalid step parameter**: Use valid step names (mapping, groupreadsbyumi, unmapped_consensus, filterconsensus, calling, allmoleculesfile)
6. **Incompatible input format**: Check column requirements for each step
7. **Missing BAM indices**: Calling step requires both BAM and CSI files
8. **Incorrect file format**: Ensure BAM files match the expected processing stage

### Error Messages:

- `VALIDATION ERROR: No VCF files generated` → Check input file paths and format
- `REFERENCE ERROR: Expected reference file missing` → Verify test data completeness
- `PRECISION FAILURE: VCF precision below threshold` → Check algorithm parameters
- `SAMPLESHEET ERROR: Invalid step choice` → Use valid step parameter values
- `INPUT ERROR: Missing required columns` → Check CSV header matches step requirements
- `FILE ERROR: BAM index not found` → Ensure CSI files exist for calling step

---

For additional support or questions about input configurations, consult the main deepUMIcaller documentation or contact the development team.
