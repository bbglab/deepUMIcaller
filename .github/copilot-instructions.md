# deepUMIcaller Pipeline Development Instructions

Always reference these instructions first and fallback to search or bash commands only when you encounter unexpected information that does not match the info here.

## Pipeline Overview

**deepUMIcaller** is a Nextflow bioinformatics pipeline for duplex sequencing data analysis that:
- Processes paired-end FASTQ files from duplex sequencing libraries
- Extracts UMI (Unique Molecular Identifier) information
- Generates duplex consensus reads for improved accuracy
- Performs variant calling with VarDict
- Provides comprehensive quality control and metrics
- Outputs VCF files and processed BAM files ready for downstream analysis

The pipeline implements the [fgbio Best Practices FASTQ to Consensus Pipeline](https://github.com/fulcrumgenomics/fgbio/blob/main/docs/best-practice-consensus-pipeline.md) with additional quality control and variant calling capabilities.

## Working Effectively

### Bootstrap and Setup
- Install Nextflow (>=25.04.1):
  ```bash
  wget -qO- https://get.nextflow.io | bash
  chmod +x nextflow
  # OR download directly:
  wget https://github.com/nextflow-io/nextflow/releases/download/v25.04.5/nextflow
  chmod +x nextflow
  ```
- Verify installation: `./nextflow -version`
- Install containerization engine (REQUIRED - choose one):
  - Docker: `apt-get install docker.io` (already available in most CI/dev environments)
  - Singularity: Follow [installation guide](https://singularity-tutorial.github.io/01-installation/)
  - Alternative: Podman, Shifter, or Charliecloud
- Install Python dependencies for utility scripts:
  ```bash
  pip install pandas numpy matplotlib seaborn pysam click
  ```

### Pipeline Validation and Testing
- **CRITICAL**: NEVER CANCEL builds or pipeline runs. Bioinformatics pipelines can take 45+ minutes for full datasets.
- Test with minimal dataset: 
  ```bash
  ./nextflow run . -profile test,docker --outdir test_results/ --max_cpus 2 --max_memory 4.GB
  ```
- **TIMING**: Test profile takes ~15-30 minutes. NEVER CANCEL. Set timeout to 60+ minutes.
- Full test dataset:
  ```bash  
  ./nextflow run . -profile test_full,docker --outdir full_test_results/
  ```
- **TIMING**: Full test takes 45+ minutes to several hours. NEVER CANCEL. Set timeout to 180+ minutes.

### Pipeline Execution Modes
The pipeline supports 4 entry points. Always specify the entry point with `--step`:

1. **Complete pipeline (default)**: `--step mapping`
   ```bash
   ./nextflow run . -profile docker \
     --input input.csv \
     --ref_fasta hs38DH.fa \
     --targetsfile targets.bed \
     --outdir results/ \
     --step mapping
   ```

2. **From UMI grouping**: `--step groupreadsbyumi`
   ```bash
   ./nextflow run . -profile docker \
     --input input.groupbyumi.csv \
     --ref_fasta hs38DH.fa \
     --targetsfile targets.bed \
     --outdir results/ \
     --step groupreadsbyumi
   ```

3. **From consensus filtering**: `--step filterconsensus`
   ```bash
   ./nextflow run . -profile docker \
     --input input.filterconsensus.csv \
     --ref_fasta hs38DH.fa \
     --targetsfile targets.bed \
     --outdir results/ \
     --step filterconsensus
   ```

4. **Variant calling only**: `--step calling`
   ```bash
   ./nextflow run . -profile docker \
     --input input.calling.csv \
     --ref_fasta hs38DH.fa \
     --targetsfile targets.bed \
     --outdir results/ \
     --step calling
   ```

### Required Input Files
- **Sample sheet**: CSV format with different columns per entry point:
  - **mapping**: `sample,fastq_1,fastq_2,read_structure`
  - **groupreadsbyumi**: `sample,bam`
  - **filterconsensus**: `sample,bam,csi,read`
  - **calling**: `sample,duplexbam,csi`
- **Reference FASTA**: Must include BWA index files in the same directory (`.bwt`, `.pac`, `.ann`, `.amb`, `.sa`)
- **Targets BED file**: Regions of interest for analysis (3-column BED format)
- See examples: `assets/samplesheet.csv`, `assets/input_example.csv`, `assets/input_example.filterconsensus_step.csv`

### Read Structure Format
Use [fgbio read structure format](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures):
- `8M1S+T 8M1S+T` - 8bp molecular barcode, 1bp skip, template to end (both reads)
- `10M1S+T 10M1S+T` - 10bp molecular barcode, 1bp skip, template to end
- `12M+T +T` - 12bp molecular barcode read 1, no barcode read 2

### Configuration and Resource Requirements
- **First run or 50+ samples**: Request minimum 8 CPUs and 20 GB memory
- **Small datasets**: 4 CPUs and 8 GB memory sufficient
- **Container pulling**: First run downloads containers (~10-15 minutes). NEVER CANCEL.
- Use appropriate profile: `test`, `docker`, `singularity`, `podman`, `conda` (conda as last resort)

### Confidence Levels
Control output detail with confidence parameters:
- `--duplex_high_conf true/false` (default: false)
- `--duplex_med_conf true/false` (default: true) 
- `--duplex_low_conf true/false` (default: false)

## Validation and Quality Control

### Always Run These Commands Before Committing
- Validate samplesheet format: 
  ```bash
  python bin/check_samplesheet.py assets/samplesheet.csv /tmp/validated.csv --step mapping
  ```
- Basic syntax check: `./nextflow config` (requires network for first-time setup)
- Test profile: `./nextflow run . -profile test,docker --outdir test_output/`
- **TIMING**: Allow 60+ minutes for test completion. NEVER CANCEL.
- **NOTE**: Nextflow requires internet access for initial setup and container pulling

### Manual Validation Scenarios
After making changes, ALWAYS run through these scenarios:

1. **Test complete pipeline**:
   ```bash
   ./nextflow run . -profile test,docker --outdir validation_test/
   ```
   - Verify output directories created: `multiqc/`, `mutations_vcf/`, `sortbamduplexconsmed/`
   - Check main outputs: `<sample>.med.vcf`, `<sample>.bam`
   - Review MultiQC report: `multiqc/multiqc_report.html`

2. **Test entry point flexibility**:
   ```bash
   ./nextflow run . -profile test,docker --step calling --outdir calling_test/
   ```
   - Verify variant calling runs independently
   - Check VCF outputs are generated

3. **Test resource scaling**:
   ```bash
   ./nextflow run . -profile test,docker --max_cpus 1 --max_memory 2.GB --outdir resource_test/
   ```
   - Ensure pipeline adapts to limited resources

4. **Test samplesheet validation**:
   ```bash
   python bin/check_samplesheet.py assets/samplesheet.csv /tmp/test_output.csv --step mapping
   python bin/check_samplesheet.py assets/input_example.filterconsensus_step.csv /tmp/test_fc.csv --step filterconsensus
   ```
   - Verify samplesheet formats are correctly validated

### Expected Output Structure
```
results/
├── multiqc/                    # QC summary report
│   └── multiqc_report.html
├── mutations_vcf/              # Variant calls
│   └── <sample>.med.vcf
├── sortbamduplexconsmed/       # Consensus BAM files
│   └── <sample>.bam
├── familymetrics/              # UMI family size metrics
├── collectduplexseqmetrics/    # Duplex sequencing QC
├── fastqc/                     # FastQ quality control
├── pipeline_info/              # Execution reports
└── [additional confidence-specific directories if enabled]
```

### Expected Build/Run Times
- **Container download**: 10-15 minutes (first time only). NEVER CANCEL.
- **Test profile**: 15-30 minutes. NEVER CANCEL. Use 60+ minute timeout.
- **Small dataset (1-5 samples)**: 30-60 minutes. NEVER CANCEL. Use 120+ minute timeout.
- **Medium dataset (10-50 samples)**: 2-4 hours. NEVER CANCEL. Use 300+ minute timeout.
- **Large dataset (50+ samples)**: 4-12 hours. NEVER CANCEL. Use 720+ minute timeout.

## Development and Debugging

### Key Code Locations
- **Main workflow**: `workflows/deepumicaller.nf`
- **Entry point**: `main.nf`
- **Module definitions**: `modules/local/` and `modules/nf-core/`
- **Configuration**: `nextflow.config` and `conf/`
- **Python utilities**: `bin/` (samplesheet validation, VCF processing, plotting)
- **Test data**: `assets/` (example inputs)

### Common Development Tasks
- **Adding new parameters**: Update `nextflow_schema.json` and `nextflow.config`
- **Module modifications**: Edit files in `modules/local/` and test with specific process
- **Configuration changes**: Modify `conf/modules.config` for process-specific settings
- **Python script changes**: Edit `bin/` scripts and test standalone

### Dependencies and Tools Used
- **Core tools**: fgbio, BWA, VarDict, SAMtools, BEDtools
- **Python dependencies**: pandas, numpy, matplotlib, seaborn, pysam, click
- **Container registries**: biocontainers (quay.io), Galaxy depot
- **QC tools**: FastQC, MultiQC, Qualimap

### Troubleshooting Common Issues
- **Out of memory**: Increase `--max_memory` parameter or reduce `--max_cpus`
- **Container issues**: Verify Docker/Singularity is running and accessible
- **Network timeouts**: Container pulls may fail; retry or use cached images
- **Missing reference files**: Ensure BWA index exists alongside FASTA file
  - Required: `.bwt`, `.pac`, `.ann`, `.amb`, `.sa` files
- **Invalid samplesheet**: Run `python bin/check_samplesheet.py input.csv output.csv --step <step>` to validate
- **Read structure errors**: Verify format matches [fgbio documentation](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)
- **Python dependency errors**: Install with `pip install pandas numpy matplotlib seaborn pysam click`
- **Nextflow download fails**: Requires internet access for initial setup and dependency downloads

## Repository Structure Reference
```
├── .github/                    # GitHub configurations
├── assets/                     # Example inputs and schemas
├── bin/                        # Python utility scripts
├── conf/                       # Configuration profiles
├── docs/                       # Documentation and images
├── modules/                    # Pipeline modules
│   ├── local/                  # Custom modules
│   └── nf-core/               # Community modules  
├── subworkflows/              # Reusable workflow components
├── workflows/                 # Main workflow definition
├── main.nf                    # Pipeline entry point
├── nextflow.config           # Main configuration
└── nextflow_schema.json      # Parameter schema
```

## Critical Reminders
- **NEVER CANCEL** any build or pipeline run operation
- Always use appropriate timeouts (60+ minutes for builds, 120+ minutes for tests)
- Test all changes with the test profile before committing
- Verify container compatibility across Docker/Singularity profiles
- Check output file generation and MultiQC reports after runs
- Document new parameters in nextflow_schema.json
- This is a bioinformatics pipeline - execution times are expected to be long

## Quick Reference Commands

### Essential Validation Commands
```bash
# Validate samplesheet
python bin/check_samplesheet.py assets/samplesheet.csv /tmp/output.csv --step mapping

# Quick test run (15-30 minutes)
./nextflow run . -profile test,docker --outdir test_results/

# Full pipeline run
./nextflow run . -profile docker \
  --input samplesheet.csv \
  --ref_fasta reference.fa \
  --targetsfile targets.bed \
  --outdir results/
```

### Common Development Tasks
```bash
# Check configuration syntax
./nextflow config

# List available profiles
./nextflow run . --help

# Resume failed run
./nextflow run . -resume -profile docker --input samplesheet.csv --outdir results/

# Clean work directory
rm -rf work/ .nextflow*
```