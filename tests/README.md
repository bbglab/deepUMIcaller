# deepUMIcaller Test Suite Documentation

## Overview

This directory contains a comprehensive test suite for the deepUMIcaller pipeline, designed to validate functionality across different input configurations and parameter combinations.

## Test Structure

```
tests/
├── main.nf.test           # Main test suite (nf-test format)
├── test_data/
│   ├── input/             # Test input CSV files and FASTQ data
│   └── expected_output/   # Reference VCF files for validation
├── nextflow.config        # Test-specific configuration
└── README.md             # This documentation
```

## Running Tests

### Prerequisites
- Nextflow installed with nf-test plugin
- deepUMIcaller dependencies available
- Test data properly configured

### Commands

```bash
# Run all tests
nf-test test

# Run specific test by tag
nf-test test --tag "normal"
nf-test test --tag "multi-file"
nf-test test --tag "split_by_chrom"

# Run with verbose output
nf-test test --verbose

# Run and update snapshots (if using snapshot testing)
nf-test test --update-snapshot
```

## Test Categories

### 1. Basic Functionality Tests

#### Test 1: Basic Functionality - Standard Processing
- **Tag**: `normal`
- **Purpose**: Validates core pipeline functionality with default parameters
- **Input**: Single FASTQ pair per sample
- **Expected**: One VCF per sample, standard processing workflow

#### Test 2: Performance Optimization - Chromosome-based Parallelization  
- **Tag**: `split_by_chrom`
- **Purpose**: Validates chromosome splitting for performance optimization
- **Input**: Same as Test 1, but with chromosome splitting enabled
- **Expected**: Same results as Test 1 but with parallel chromosome processing

### 2. Multi-file Handling Tests

#### Test 3: Technical Replicates - Multi-file Sample Processing
- **Tag**: `multi-file`
- **Purpose**: Validates merging of technical replicates (sequencing lanes)
- **Input**: Multiple FASTQ pairs from the same biological sample
- **Parameters**: `splitted_original_sample = true`
- **Expected**: Merged results at sample level, improved sensitivity

#### Test 4: Biological Replicates - Multi-sample Patient Processing
- **Tag**: `multi-sample`
- **Purpose**: Validates patient-level aggregation across biological samples
- **Input**: Multiple samples with shared `parent_dna` identifier
- **Expected**: Both sample-level and patient-level VCF outputs

#### Test 5: Comprehensive Scenario - Complete Multi-level Processing
- **Tag**: `multiAll`
- **Purpose**: Most complex scenario combining all processing modes
- **Input**: Technical + biological replicates with chromosome splitting
- **Parameters**: `splitted_original_sample = true`, `split_by_chrom = true`
- **Expected**: Multi-level processing with maximum computational efficiency

## Validation Criteria

### Success Metrics
- **Pipeline Execution**: All tests must complete without errors
- **VCF Precision**: ≥99% precision when compared to reference outputs
- **File Generation**: Expected output files must be created in correct locations

### Validation Process
1. Pipeline execution validation (`workflow.success`)
2. Output file existence validation
3. VCF precision comparison using `compare_vcfs_detailed.py`
4. Detailed error reporting with specific failure reasons

## Test Data

### Input Files Location: `tests/test_data/input/`
- `input_test.csv` - Basic single-sample input
- `input_test_multi-file.csv` - Multi-lane technical replicates
- `input_test_multi-sample.csv` - Multi-sample patient data
- `input_test_multiAll.csv` - Complex combined scenario

### Reference Files Location: `tests/test_data/expected_output/`
Contains validated VCF outputs for precision comparison.

## Parameter Combinations Tested

| Test | splitted_original_sample | split_by_chrom | parent_dna | Use Case |
|------|-------------------------|----------------|------------|----------|
| 1    | false (default)         | false (default)| default    | Standard processing |
| 2    | false (default)         | true           | default    | Performance optimization |
| 3    | true                    | false (default)| default    | Technical replicates |
| 4    | false (default)         | false (default)| explicit   | Biological replicates |
| 5    | true                    | true           | explicit   | Comprehensive scenario |

## Adding New Tests

### Steps to Add a New Test:
1. **Create test data**: Add input CSV and FASTQ files to `test_data/input/`
2. **Generate reference**: Run pipeline manually and validate output
3. **Add reference data**: Place validated VCF files in `test_data/expected_output/`
4. **Write test case**: Add new test block to `main.nf.test`
5. **Document test**: Update this README with test description
6. **Validate**: Run new test to ensure it passes

### Test Documentation Template:
```nextflow
/*
========================================================================================
    TEST X: CATEGORY - Descriptive Name
========================================================================================
Purpose: Clear description of what this test validates

Input: Description of test data and structure
Parameters: Key parameter settings and their purpose

Expected Behavior:
    - Step 1: What should happen first
    - Step 2: What should happen next
    - Step N: Final expected outcome

Use Case: When this configuration would be used in practice

Success Criteria: Specific metrics for success
========================================================================================
*/
test("Descriptive test name") {
    tag "category"
    
    when {
        params {
            // Commented parameter explanations
        }
    }
    
    then {
        // Clear assertions with descriptive error messages
    }
}
```

## Troubleshooting

### Common Issues:
1. **Missing test data**: Ensure FASTQ files are properly linked/available
2. **Reference mismatches**: Update expected outputs when pipeline changes
3. **Precision failures**: Investigate algorithm changes or parameter effects
4. **Resource limitations**: Some tests may require significant compute resources

### Debug Tips:
- Use `--verbose` flag for detailed output
- Check intermediate files in work directories
- Compare outputs manually using `compare_vcfs_detailed.py`
- Validate test data integrity before running tests

## Contributing

When modifying the pipeline:
1. **Run full test suite** before submitting changes
2. **Update tests** if functionality changes
3. **Add new tests** for new features
4. **Update documentation** to reflect changes

## Contact

For test-related questions or issues, consult the main deepUMIcaller documentation or contact the development team.
