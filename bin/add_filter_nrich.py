#!/usr/bin/env python

# usage: python3 add_filter_nrich.py vcf_file ns_position_file output_filename filter_name
# use filter_name = "n_rich"
# ns_position_file generated from the extended exons file to avoid overdistributed values; also valid to use the panel coordinates

import pandas as pd
import numpy as np
import sys

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


def add_filter(old_filt, add_filt, filt_name):
        """
        old filt is the current FILTER field value
        add_filt is a boolean, either True or False
        filt_name is the name that should be added in the FILTER field in case the add_filt value is True
        """
        if add_filt:
            if old_filt == "PASS":
                return filt_name
            old_filt += ";" + filt_name
        return ";".join( sorted(old_filt.split(";")) )


def add_filter_nrich(vcf, ns_position_file, filter_name, min_depth_valid = 5):
    """
    Adds to a VCF-like dataframe an additional filter
    in mutations located at positions which have more 
    Ns than expected. The threshold is defined as 
    the median+2*s.d.
    
    Parameters
    ----------
    vcf: pandas DataFrame
        VCF-like dataframe containing the mutations

    Returns
    -------
    updated_vcf: pandas DataFrame
        VCF-like dataframe updated
    """

    # define filter threshold
    ## load file with Ns counts per position along the entire panel (no mutations only)
    ns_position_df = pd.read_csv(ns_position_file, sep = "\t", header = None,
                                names = ["CHROM", "START", "TOTAL_DEPTH", "N_COUNT"])
    print(ns_position_df.shape)

    over_valid_depth = ns_position_df["TOTAL_DEPTH"][ns_position_df["TOTAL_DEPTH"] > min_depth_valid]
    print(over_valid_depth.shape)
    
    # TODO: add more documentation on this
    # we use positions covered by at least half of the median coverage
    min_depth_threshold = over_valid_depth.quantile(0.75) #* 0.5
    print(min_depth_threshold)
    del over_valid_depth

    ns_position_df = ns_position_df[ns_position_df["TOTAL_DEPTH"] > min_depth_threshold].reset_index(drop = True)
    print(ns_position_df.shape)

    ## calculate the proportion of Ns per position adding a pseudocount to enable log-transformation
    ns_position_df["N_COUNT/TOTAL_DEPTH_pseudocount"] = ns_position_df.apply(lambda row: (row["N_COUNT"] / row["TOTAL_DEPTH"])+0.0000001,
                                                                                axis = 1)
    ## compute threshold over the log-transformed data and pass back to the original scale
    log_data = np.log(ns_position_df["N_COUNT/TOTAL_DEPTH_pseudocount"])
    median = np.median(log_data)
    std_dev = np.std(log_data)
    threshold = np.exp(median+2*std_dev)
    
    
    # add filter when the proportion of Ns in the mutated position is higher than the threshold
    #(w/ pseudocount to fit with the threshold calculation)
        #NDP: 9th column in last field of VCF (0-based)
        #CDP: 7th column in last field of VCF (0-based)
    def annot(sample, filter, threshold):

        # if there is no_pileup_support, the filter is not annotated because depths are zero
        if "no_pileup_support" not in filter:

            proportion_ns = (int(sample.split(":")[9]) / (int(sample.split(":")[9])+int(sample.split(":")[7])))+0.0000001

            if proportion_ns > threshold:
                return True

        return False

    # compute which mutations have a proportion of ns higher than the threshold
    vcf[filter_name] = vcf.apply(lambda row: annot(row["SAMPLE"], row["FILTER"], threshold), axis = 1)

    # add new value to FILTER column when needed
    vcf["FILTER"] = vcf[["FILTER",filter_name]].apply(
                                                        lambda x: add_filter(x["FILTER"], x[filter_name], filter_name),
                                                        axis = 1
                                                        )

    # remove unneeded column from vcf
    return vcf.drop([filter_name], axis = 1), threshold


def main(vcf_file, ns_position_file, output_filename, filter_name, min_valid_depth):
    """
    Reads a VCF file and adds n_rich value to the
    FILTER field where applicable.
    
    Parameters
    ----------
    vcf_file: str:
        Path to the VCF file to be updated.
    ns_position_file: str
        Path to the file containing the number of 
        Ns per position.
    output_filename: str
        Name path for the resulting updated VCF.
    """

    ###
    # Read the VCF file body and add 
    ###
    vcf = pd.read_csv(vcf_file, sep = '\t', header = None, comment= '#')
    vcf.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]

    ###
    # Add filter: mutations in positions with more Ns than expected
    ###
    updated_vcf, threshold = add_filter_nrich(vcf, ns_position_file, filter_name, min_valid_depth)

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
    header_lines.append(f'##FILTER=<ID={filter_name},Description="Variant located in a position with more Ns than expected. Threshold in this sample: {round(threshold, 6)}. deepUMIcaller.">')
    
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
    ns_position_file = sys.argv[2]
    output_filename = sys.argv[3]
    filter_name = sys.argv[4]
    min_valid_depth = int(sys.argv[5])
    main(vcf_file, ns_position_file, output_filename, filter_name, min_valid_depth)