"""
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCF Comparison Tool for deepUMIcaller Pipeline Testing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Description:
    This script compares two VCF (Variant Call Format) files and calculates various
    similarity metrics. It's primarily used in the deepUMIcaller pipeline testing
    framework to validate that pipeline outputs match expected reference results.

Purpose:
    - Quality assurance for variant calling pipelines
    - Regression testing to ensure pipeline changes don't affect results
    - Performance evaluation using multiple similarity metrics
    - Debugging variant calling differences

Metrics Calculated:
    1. Jaccard Index: |intersection| / |union| * 100
       - Measures overall similarity between two variant sets
       - Range: 0-100% (100% = identical variant sets)
    
    2. Precision: |intersection| / |vcf1_variants| * 100
       - Percentage of VCF1 variants that are also in VCF2
       - Answers: "How many of my results are correct?"
    
    3. Recall: |intersection| / |vcf2_variants| * 100
       - Percentage of VCF2 variants that are also in VCF1
       - Answers: "How many of the expected results did I find?"
    
    4. F1 Score: 2 * (precision * recall) / (precision + recall)
       - Harmonic mean of precision and recall
       - Balanced measure when you care about both false positives and false negatives

Usage Examples:
    # Basic comparison (returns Jaccard index)
    python compare_vcfs_detailed.py current.vcf expected.vcf
    
    # Get precision score (used in nf-test validation)
    python compare_vcfs_detailed.py current.vcf expected.vcf precision
    
    # Debug mode - detailed comparison with variant lists
    python compare_vcfs_detailed.py current.vcf expected.vcf debug

Integration with nf-test:
    This script is called from the test validation functions in main.nf.test:
    
    def cmd = "python3 bin/compare_vcfs_detailed.py ${currentVcf} ${expectedVcf}"
    def precision = proc.in.text.trim() as Double
    assert precision >= 99.0 : "PRECISION FAILURE: VCF precision below threshold"

Input Requirements:
    - VCF files must be properly formatted with standard columns
    - Minimum required columns: CHROM, POS, ID, REF, ALT
    - Empty lines and header lines (starting with #) are automatically skipped
    - Malformed lines are silently skipped to ensure robust parsing

Output Formats:
    - Default modes: Returns integer percentage (0-100)
    - Debug mode: Returns detailed comparison statistics and variant lists
    - Designed for both human reading (debug) and automated parsing (default)

Author: deepUMIcaller pipeline team
Version: Compatible with VCF format specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""

import sys
import re

def load_variants(vcf_file):
    """
    Load variants from a VCF file and return them as a set of tuples.
    
    Args:
        vcf_file (str): Path to the VCF file to parse
        
    Returns:
        tuple: (variants_set, variant_details_dict)
            - variants_set: Set of (chrom, pos, ref, alt) tuples
            - variant_details_dict: Dictionary mapping variant keys to full VCF lines
            
    Description:
        Parses a VCF file and extracts variant information. Each variant is represented
        as a tuple of (chromosome, position, reference, alternate) to enable easy
        comparison between different VCF files. Skips header lines, empty lines,
        and malformed entries to ensure robust parsing.
        
    VCF Format Expected:
        Standard VCF format with tab-separated columns:
        CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  [SAMPLES...]
        
    Error Handling:
        - Skips empty lines gracefully
        - Skips header lines (starting with #)
        - Skips malformed lines with insufficient columns
        - No exceptions thrown for parsing errors
    """
    variants = set()
    variant_details = {}
    with open(vcf_file, "r") as f:
        for line in f:
            if not line.strip():  # skip empty lines
                continue
            if line.startswith("#"):  # skip header
                continue

            cols = re.split(r"\s+", line.strip())
            if len(cols) < 5:  # malformed line, skip
                continue

            chrom = cols[0]
            pos = cols[1]
            ref = cols[3]
            alt = cols[4]
            variant_key = (chrom, pos, ref, alt)
            variants.add(variant_key)
            variant_details[variant_key] = line.strip()
    return variants, variant_details

def main(vcf1, vcf2, mode="jaccard"):
    """
    Compare two VCF files and calculate similarity metrics.
    
    Args:
        vcf1 (str): Path to the first VCF file (typically current/test results)
        vcf2 (str): Path to the second VCF file (typically expected/reference results)
        mode (str): Comparison mode - "jaccard", "precision", "recall", "f1", or "debug"
        
    Description:
        Performs comprehensive comparison between two VCF files using set theory.
        The comparison is based on exact matching of (chromosome, position, reference, alternate)
        tuples, ensuring that variants are identical in all key fields.
        
    Modes Explained:
        - jaccard: Overall similarity measure (default for general comparison)
        - precision: Used in pipeline testing to ensure accuracy of results
        - recall: Measures completeness of variant detection
        - f1: Balanced measure for comprehensive evaluation
        - debug: Detailed output for troubleshooting and analysis
        
    Mathematical Formulations:
        Let A = variants in vcf1, B = variants in vcf2
        - Jaccard = |A ∩ B| / |A ∪ B| × 100
        - Precision = |A ∩ B| / |A| × 100
        - Recall = |A ∩ B| / |B| × 100
        - F1 = 2 × (Precision × Recall) / (Precision + Recall)
        
    Pipeline Integration:
        In deepUMIcaller testing, this function is called with mode="precision"
        to validate that test results achieve ≥99% precision against reference data.
        This ensures that pipeline changes don't introduce false positive variants.
        
    Error Conditions:
        - Prints "ERROR: No variants in vcf1" if first VCF is empty
        - Returns 0 for all metrics if no variants found
        - Handles division by zero gracefully
    """
    v1_variants, v1_details = load_variants(vcf1)
    v2_variants, v2_details = load_variants(vcf2)

    # Validate input data
    if not v1_variants:
        print("ERROR: No variants in vcf1")
        return

    # Calculate set operations for comparison metrics
    overlap = v1_variants.intersection(v2_variants)  # Variants present in both files
    union = v1_variants.union(v2_variants)          # All unique variants across both files
    
    # Calculate different similarity metrics
    jaccard = (len(overlap) / len(union)) * 100 if union else 0        # Overall similarity
    precision = (len(overlap) / len(v1_variants)) * 100 if v1_variants else 0  # Accuracy of vcf1
    recall = (len(overlap) / len(v2_variants)) * 100 if v2_variants else 0     # Completeness vs vcf2
    f1_score = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0  # Balanced measure

    # Output based on requested mode
    if mode == "debug":
        # Comprehensive debugging information for analysis and troubleshooting
        print(f"VCF1 variants: {len(v1_variants)}")
        print(f"VCF2 variants: {len(v2_variants)}")
        print(f"Overlap: {len(overlap)}")
        print(f"Union: {len(union)}")
        print(f"Jaccard Index: {jaccard:.2f}%")
        print(f"Precision (vcf1 in vcf2): {precision:.2f}%")
        print(f"Recall (vcf2 in vcf1): {recall:.2f}%")
        print(f"F1 Score: {f1_score:.2f}%")
        
        # Show variants only in vcf1 (false positives if vcf2 is ground truth)
        only_v1 = v1_variants - v2_variants
        if only_v1:
            print(f"\nVariants only in VCF1 ({len(only_v1)}) - Potential False Positives:")
            for var in sorted(list(only_v1))[:10]:  # Show first 10 for readability
                print(f"  {var[0]}:{var[1]} {var[2]}>{var[3]}")
            if len(only_v1) > 10:
                print(f"  ... and {len(only_v1) - 10} more variants")
        
        # Show variants only in vcf2 (false negatives if vcf2 is ground truth)
        only_v2 = v2_variants - v1_variants
        if only_v2:
            print(f"\nVariants only in VCF2 ({len(only_v2)}) - Potential False Negatives:")
            for var in sorted(list(only_v2))[:10]:  # Show first 10 for readability
                print(f"  {var[0]}:{var[1]} {var[2]}>{var[3]}")
            if len(only_v2) > 10:
                print(f"  ... and {len(only_v2) - 10} more variants")
                
    elif mode == "precision":
        # Used in automated testing - returns integer percentage for threshold comparison
        print(int(precision))
    elif mode == "recall":
        # Alternative metric for completeness evaluation
        print(int(recall))
    elif mode == "f1":
        # Balanced metric combining precision and recall
        print(int(f1_score))
    else:  # default jaccard mode
        # General similarity measure - useful for overall comparison
        print(int(jaccard))

if __name__ == "__main__":
    """
    Command-line interface for the VCF comparison tool.
    
    Usage:
        python compare_vcfs_detailed.py <vcf1> <vcf2> [mode]
        
    Arguments:
        vcf1: Path to first VCF file (current/test results)
        vcf2: Path to second VCF file (expected/reference results)
        mode: Optional comparison mode (default: jaccard)
              Valid modes: jaccard, precision, recall, f1, debug
              
    Examples:
        # Basic comparison
        python compare_vcfs_detailed.py test.vcf reference.vcf
        
        # Pipeline testing (used in nf-test)
        python compare_vcfs_detailed.py current.vcf expected.vcf precision
        
        # Detailed analysis
        python compare_vcfs_detailed.py test.vcf reference.vcf debug
        
    Exit Codes:
        0: Successful execution
        1: Invalid arguments or usage error
    """
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print(f"Usage: python {sys.argv[0]} <vcf1> <vcf2> [mode]")
        print("Modes: jaccard (default), precision, recall, f1, debug")
        sys.exit(1)

    mode = sys.argv[3] if len(sys.argv) == 4 else "jaccard"
    main(sys.argv[1], sys.argv[2], mode)
