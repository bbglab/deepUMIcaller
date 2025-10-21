#!/usr/bin/env python3

import sys
import argparse
from collections import defaultdict

def validate_vcf_biology(vcf_file):
    """
    Validate that VCF contains biologically reasonable variants
    """
    issues = []
    variant_types = defaultdict(int)
    allele_frequencies = []
    
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            if not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 8:
                continue
                
            chrom, pos, id, ref, alt, qual, filter_status, info = fields[:8]
            
            # Check chromosome names are valid
            if not chrom.startswith('chr') and not chrom.isdigit():
                issues.append(f"Invalid chromosome name: {chrom}")
            
            # Check position is numeric
            try:
                pos_int = int(pos)
                if pos_int <= 0:
                    issues.append(f"Invalid position: {pos}")
            except ValueError:
                issues.append(f"Non-numeric position: {pos}")
            
            # Check REF and ALT are valid DNA bases
            valid_bases = set('ATCGN')
            if not set(ref.upper()).issubset(valid_bases):
                issues.append(f"Invalid REF allele: {ref}")
            if not set(alt.upper()).issubset(valid_bases):
                issues.append(f"Invalid ALT allele: {alt}")
            
            # Classify variant type
            if len(ref) == 1 and len(alt) == 1:
                variant_types['SNV'] += 1
            elif len(ref) > len(alt):
                variant_types['Deletion'] += 1
            elif len(ref) < len(alt):
                variant_types['Insertion'] += 1
            else:
                variant_types['Complex'] += 1
            
            # Extract allele frequency if available
            if 'AF=' in info:
                try:
                    af = float(info.split('AF=')[1].split(';')[0])
                    allele_frequencies.append(af)
                except:
                    pass
    
    # Biological reasonableness checks
    total_variants = sum(variant_types.values())
    if total_variants == 0:
        issues.append("No variants found - this might indicate a problem")
    
    # Check variant type distribution
    snv_ratio = variant_types['SNV'] / total_variants if total_variants > 0 else 0
    if snv_ratio < 0.5:
        issues.append(f"Low SNV ratio ({snv_ratio:.2f}) - expected majority to be SNVs")
    
    # Check allele frequencies are reasonable
    if allele_frequencies:
        high_af_count = sum(1 for af in allele_frequencies if af > 0.9)
        if high_af_count / len(allele_frequencies) > 0.1:
            issues.append("Too many variants with very high allele frequency (>90%)")
    
    return {
        'valid': len(issues) == 0,
        'issues': issues,
        'stats': {
            'total_variants': total_variants,
            'variant_types': dict(variant_types),
            'avg_allele_freq': sum(allele_frequencies) / len(allele_frequencies) if allele_frequencies else 0
        }
    }

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate VCF biological correctness')
    parser.add_argument('vcf', help='VCF file to validate')
    parser.add_argument('--verbose', action='store_true', help='Show detailed output')
    
    args = parser.parse_args()
    
    result = validate_vcf_biology(args.vcf)
    
    if args.verbose:
        print(f"Total variants: {result['stats']['total_variants']}")
        print(f"Variant types: {result['stats']['variant_types']}")
        print(f"Average AF: {result['stats']['avg_allele_freq']:.3f}")
        
        if result['issues']:
            print("\nIssues found:")
            for issue in result['issues']:
                print(f"  - {issue}")
    
    # Exit code: 0 if valid, 1 if issues found
    sys.exit(0 if result['valid'] else 1)
