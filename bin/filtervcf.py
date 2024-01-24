#!/usr/bin/env python


import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

######
# Read VCF file coming from VarDict2
#      Note that the file can only contain one sample
#      if two samples are given in the same VCF it will remove the first sample in the file
######

def read_from_vardict_VCF_all(sample,
                                name,
                                subset_val = 0.35, 
                                columns_to_keep = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE",
                                                    'DEPTH', 'REF_DEPTH', 'ALT_DEPTH', 'VAF',
                                                    'vd_DEPTH', 'vd_REF_DEPTH', 'vd_ALT_DEPTH'],
                                n_bins = 100,
                                location = "",
                                plottingDist = True,
                                keep_original_VAF = False
                                ):
    """
    Read VCF file coming from Vardict
    Note that the file can only contain one sample
    if two samples are given in the same VCF it will remove the first sample in the file
    
    Mandatory arguments:
        sample,
        name,
    
    Optional arguments:
        subset_val = 0.35, 
        columns_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'ALT_DEPTH', 'VAF'], # add 'PID' for phased mutations
        n_bins = 100,
        location = "",
        plottingDist = True,
        only_SNVs = True
    """

    dat = pd.read_csv(location + name,
                    sep = '\t', header = None, comment= '#')
    dat.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]


    formats = dat.loc[:,"FORMAT"].str.split(':')
    values = dat.loc[:,"SAMPLE"].str.split(':')

    features_list = list()
    for lab, val in zip(formats, values):
        features_list.append({ l:v for l, v in zip(lab, val)})
    features_df = pd.DataFrame(features_list)

    dat_full = pd.concat((dat, features_df), axis = 1)

    if "AD" not in dat_full.columns:
        print(f"AD not present in the format field, revise the VCF reading function")
        return dat_full
    else:
        
        for ele in ["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]:
            if ele not in dat_full.columns:
                print(f"{ele} not present in the format field, revise the VCF reading function")
                return dat_full
        
        ## variant depth computation
        variant_depth_precomputed = dat_full["VD"].astype(int)
        variant_depth_from_splitted = [int(v[1]) for v in dat_full["AD"].str.split(",")]
        variant_depth_from_strands = [int(v[0]) + int(v[1]) for v in dat_full["ALD"].str.split(",")]
        
        variant_depth_max1 = [ max(dp1, dp2) for dp1, dp2 in zip(variant_depth_precomputed, variant_depth_from_splitted) ]
        variant_depth_max2 = [ max(dp1, dp2) for dp1, dp2 in zip(variant_depth_max1, variant_depth_from_strands) ]
        
        # assign it to the column
        dat_full["vd_ALT_DEPTH"] = variant_depth_max2
        
        ## reference depth computation
        reference_depth_from_splitted = [int(v[0]) for v in dat_full["AD"].str.split(",")]
        reference_depth_from_strands = [int(v[0]) + int(v[1]) for v in dat_full["RD"].str.split(",")]
        reference_depth_max = [ max(dp1, dp2) for dp1, dp2 in zip(reference_depth_from_splitted, reference_depth_from_strands) ]
        
        dat_full["vd_REF_DEPTH"] = reference_depth_max
        
        
        ## add the maximum of one side and the other
        total_depth_from_both = [int(dp1) + int(dp2) for dp1, dp2 in zip(reference_depth_max, variant_depth_max2) ]
        
        ## total depth computation from the data in the vcf
        total_depth_precomputed = dat_full["DP"].astype(int)
        total_depth_from_splitted = [int(v[0]) + int(v[1]) for v in dat_full["AD"].str.split(",")]
        
        total_depth_max1 = [ max(dp1, dp2) for dp1, dp2 in zip(total_depth_precomputed, total_depth_from_splitted) ]
        
        ## compute the maximum total depth from the two options computed
        total_depth_max2 = [ max(dp1, dp2) for dp1, dp2 in zip(total_depth_from_both, total_depth_max1) ]
        
        # assign the total depth to the DEPTH column
        dat_full["vd_DEPTH"] = total_depth_max2



    ##
    # Corrected depths
    ##
    for ele in ["CDP", "CAD", "NDP"]:
        if ele not in dat_full.columns:
            print(f"{ele} not present in the format field, revise the VCF reading function")
            return dat_full

    # assign it to the column
    dat_full["ALT_DEPTH"] = [int(v[1]) for v in dat_full["CAD"].str.split(",")]
    dat_full["REF_DEPTH"] = [int(v[0]) for v in dat_full["CAD"].str.split(",")]

    dat_full["DEPTH"] = dat_full["CDP"].astype(int)

    dat_full["numNs"] = dat_full["NDP"].astype(int)


    if "AF" not in formats or not keep_original_VAF:
        dat_full["VAF"] = dat_full["ALT_DEPTH"] / dat_full["DEPTH"]
    
    else:
        dat_full["VAF"] = dat_full["AF"].astype(float)
    
    
    if plottingDist:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = [14, 4])
        # fig, (ax1, ax2) = plt.subplots(2, 1, figsize = [13, 16])
        fig.suptitle(f"Sample {sample}")

        sns.histplot(dat_full, x="VAF", bins = n_bins,
                        ax = ax1
                    )
        sns.histplot(dat_full[dat_full["VAF"] < subset_val], x="VAF", bins= n_bins,
                        ax = ax2
                    )

        # ax1.axvline(0, color = 'r', linestyle = '--')
        # ax1.axvline(subset_val, color = 'r', linestyle = '--')

        ax1.set_xlabel("VAF")
        ax2.set_xlabel("VAF")
        fig.savefig(f"{sample}.VAF_distribution.png")

    return dat_full[ columns_to_keep ].reset_index(drop = True)




def main(sample, vcf_file, filter_str, vaf_germline, output_filename):
    """
    Your script's main function.
    """

    vcf_read = read_from_vardict_VCF_all(sample, vcf_file)
    filter_list = sorted(filter_str.split(","))
    for filt in filter_list:
        print(f'{sum(vcf_read["FILTER"].str.contains(filt))} variants filtered out due to {filt} filter.')
        vcf_read = vcf_read[~vcf_read["FILTER"].str.contains(filt)]
    vcf_read = vcf_read.reset_index(drop = True)

    print(f'{sum(vcf_read["VAF"] > vaf_germline)} variants filtered for germline variants.')
    vcf_read = vcf_read[vcf_read["VAF"] <= vaf_germline].reset_index(drop = True)

    vcf_read[["CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO",
                "FORMAT", "SAMPLE"]].to_csv(output_filename, sep='\t', header=["#CHROM", "POS", "ID",
                                                                                "REF", "ALT", "QUAL",
                                                                                "FILTER", "INFO",
                                                                                "FORMAT", "SAMPLE"], index = False)

    # Print a success message or return a result if needed
    print(f"{output_filename} VCF file created successfully.")



if __name__ == '__main__':
    sample = sys.argv[1]
    vcf_file = sys.argv[2]
    filters = sys.argv[3]
    vaf_threshold = float(sys.argv[4])
    vcf_file_out = sys.argv[5]
    main(sample, vcf_file, filters, vaf_threshold, vcf_file_out)

