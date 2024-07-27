#!/usr/bin/env python

# to be added in the pipeline before executing this script to generate the vcf_lowcmappability_file
# unmappable region file: unmappable.bed
# vcf_derived_bed is the file generated from the coordinates of the mutations, extending them +5 bp in the edges and with the variant id in the format chr;pos;ref;alt
# bedtools intersect -a vcf_derived_bed -b unmappable_file -u  > unmappable_muts_file.bed

# usage: python3 add_filter_lowmappability.py vcf_file unmappable_muts_file.bed output_filename filter_name
# use filter_name = "low_mappability"

import pandas as pd
import sys
import numpy as np

def update_filter_field(row, filter_name):
    """
    Updates the FILTER columns of a VCF-like
    dataframe if filter_name has a value. 
    Respects alphabetical ordering of the values
    in the updated FILTER column
    """

    
    if row[filter_name] == filter_name:

        if row["FILTER"] != "PASS":

            filters = row["FILTER"].split(";")
            filters.append(row[filter_name])
            filters.sort(reverse = False)
            row["FILTER"] = ";".join(filters)
        
        else:
            row["FILTER"] = row[filter_name]
    
    return row

def add_filter_lowmappability(vcf, bed_low_mappability_file, filter_name):
    """
    Adds to a VCF-like dataframe an additional filter
    in mutations located at low complex or repetitive regions
    of the genome
    
    Parameters
    ----------
    vcf: pandas DataFrame
        VCF-like dataframe containing the mutations
    bed_low_mappability_file: str
        Path to the BED4 file containing the bedtools intersection
        between the VCF mutations and the low mappability regions.
        Mandatory columns: CHROM, START, END, VARIANT_ID 

    Returns
    -------
    updated_vcf: pandas DataFrame
        VCF-like dataframe updated
    """

    bed_mappability_df = pd.read_csv(bed_low_mappability_file, sep = '\t', header = None,
                                        na_filter=False)
    bed_mappability_df.columns = ["CHROM", "START", "END", "VARIANT_ID"]

    # add filter_name label for later merging with the vcf
    bed_mappability_df[filter_name] = filter_name  

    # add VARIANT_ID field to the VCF for merging
    vcf["VARIANT_ID"] = vcf.apply(lambda row: str(row["CHROM"])+";"+str(row["POS"])+";"+row["REF"]+";"+row["ALT"],
                                axis = 1)      
    vcf = vcf.merge(bed_mappability_df[["VARIANT_ID", filter_name]], how = "left", on = "VARIANT_ID")

    # add new value to FILTER column when needed
    updated_vcf = vcf.apply(lambda row: update_filter_field(row, filter_name), axis = 1)
    # remove unneeded column from vcf
    updated_vcf = updated_vcf.drop([filter_name, "VARIANT_ID"], axis = 1)

    return updated_vcf

def main(vcf_file, bed_low_mappability_file, output_filename, filter_name):
    """
    Reads a VCF file and adds low_mappability value to the
    FILTER field where applicable.
    
    Parameters
    ----------
    vcf_file: str:
        Path to the VCF file to be updated.
    bed_low_mappability_file: str
        Path to the BED4 file containing the bedtools intersection
        between the VCF mutations and the low mappability regions.
        Mandatory columns: CHROM, START, END, VARIANT_ID 
    output_filename: str
        Name path for the resulting updated VCF.
    filter_name: str
        Name to be added in the FILTER column. 
        Default: low_mappability

    """

    ###
    # Read the VCF file body and add 
    ###
    vcf = pd.read_csv(vcf_file, sep = '\t', header = None, comment= '#',
                            dtype={0: str, 1: int, 2: str, 3: str, 4: str, 5: int, 6: str, 7: str, 8: str, 9: str},
                            na_filter=False)
    vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

    ###
    # Add filter: mutations in low mappable regions
    ###
    updated_vcf = add_filter_lowmappability(vcf, bed_low_mappability_file, filter_name)


    ###
    # Read the VCF header
    ###
    first_header_lines = []
    header_lines = []
    with open(vcf_file, 'r') as vcf_open_file:
        for line in vcf_open_file:
            if line.startswith('##'):
                if line.strip().strip("#").split("=")[0].islower():
                    first_header_lines.append(line.strip())
                else:
                    header_lines.append(line.strip())
            elif line.startswith('#'):
                single_header = line.strip()
            else:
                break
    
    ###
    # Add your custom header lines
    ###    
    # FILTER field
    header_lines.append(f'##FILTER=<ID={filter_name},Description="Variant located in a genomic region with low mappability. deepUMIcaller.">')
    
    ###
    # Combine the modified header rows
    ###
    # Convert the list of header lines to a single string
    header_str = '\n'.join(sorted(first_header_lines) + sorted(header_lines))
    header_str = header_str + "\n" + single_header

    ###
    # Write the modified VCF header and the updated VCF body into a new file with the appropriate name
    ###
    with open(output_filename, 'w') as new_vcf_file:
        new_vcf_file.write(header_str + '\n')  # Write the modified header
        updated_vcf.to_csv(new_vcf_file, sep='\t', header=False, index=False)  # Write the data rows
    
    # Print a success message or return a result if needed
    print(f"{output_filename} VCF file created successfully.")

if __name__ == '__main__':
    vcf_file = sys.argv[1]
    bed_low_mappability_file = sys.argv[2]
    output_filename = sys.argv[3]
    filter_name = sys.argv[4]
    main(vcf_file, bed_low_mappability_file, output_filename, filter_name)