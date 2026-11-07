#!/usr/bin/env python


import click
import pandas as pd
import numpy as np

######
# Read VCF file coming from VarDict2
#      Note that the file can only contain one sample
#      if two samples are given in the same VCF it will remove the first sample in the file
# also it takes advantage of the fields in the FORMAT field that are coming from the recomputed depths the result from the
######


def read_from_vardict_VCF_all(sample,
                                name,
                                columns_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'REF_DEPTH', 'ALT_DEPTH', 'VAF',
                                                    'vd_DEPTH', 'vd_REF_DEPTH', 'vd_ALT_DEPTH'],
                                vaf_distortion_threshold = 3):
    """
    Read VCF file coming from Vardict2
    Note that the file can only contain one sample
    if two samples are given in the same VCF it will remove the first sample in the file
    also it takes advantage of the fields in the FORMAT field that are coming from the
    recomputed depths scripts in deepUMIcaller


    Mandatory arguments:
        sample,
        name,

    Optional arguments:
        columns_to_keep = ['CHROM', 'POS', 'REF', 'ALT', 'DEPTH', 'ALT_DEPTH', 'VAF'], # add 'PID' for phased mutations
        vaf_distortion_threshold = 3
    """

    print(f"Processing {sample}")

    dat = pd.read_csv(name,
                    sep = '\t', header = None, comment= '#', na_filter = False)
    dat.columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
    print("Total mutations:", dat.shape[0])

    formats = dat.loc[:,"FORMAT"].str.split(':')
    values = dat.loc[:,"SAMPLE"].str.split(':')

    features_list = list()
    for lab, val in zip(formats, values):
        features_list.append({ l:v for l, v in zip(lab, val)})
    features_df = pd.DataFrame(features_list)

    dat_full = pd.concat((dat, features_df), axis = 1)


    ## Note ##
    # The fields DEPTH and ALT_DEPTH might be computed differently
    #     in the current version of the code we are computing everything from the AD field.
    # dat_full["ALT_DEPTH"] = [int(v[1]) for v in dat_full["AD"].str.split(",")]
    # dat_full["DEPTH"] = dat_full["DP"].astype(int)

    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    ##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
    ##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">



    if "AD" not in dat_full.columns:
        print(f"AD not present in the format field, revise the VCF reading function")
        return dat_full

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

    # compute VAF from VarDict counts
    dat_full["vd_VAF"] = dat_full["vd_ALT_DEPTH"] / dat_full["vd_DEPTH"]





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
    dat_full["VAF_Ns"] = dat_full["numNs"] / (dat_full["DEPTH"] + dat_full["numNs"])

    # compute VAF
    dat_full["VAF"] = dat_full["ALT_DEPTH"] / dat_full["DEPTH"]

    ##
    # All molecule depths
    ##
    for ele in ["CDPAM", "CADAM", "NDPAM"]:
        if ele not in dat_full.columns:
            print(f"{ele} not present in the format field, revise your VCF")
            # return dat_full
            print(f"In the meantime, we are using {ele[:3]} to fill the information for {ele}")
            dat_full[ele] = dat_full[ele[:3]].astype(int)

    # assign it to the column
    dat_full["ALT_DEPTH_AM"] = [int(v[1]) for v in dat_full["CADAM"].str.split(",")]
    dat_full["REF_DEPTH_AM"] = [int(v[0]) for v in dat_full["CADAM"].str.split(",")]
    dat_full["DEPTH_AM"] = dat_full["CDPAM"].astype(int)

    # compute VAF_AM
    dat_full["VAF_AM"] = dat_full["ALT_DEPTH_AM"] / dat_full["DEPTH_AM"]

    dat_full["numNs_AM"] = dat_full["NDPAM"].astype(int)
    dat_full["VAF_Ns_AM"] = dat_full["numNs_AM"] / (dat_full["DEPTH_AM"] + dat_full["numNs_AM"])

    # Compute ND values
    dat_full["ALT_DEPTH_ND"] = dat_full["ALT_DEPTH_AM"] - dat_full["ALT_DEPTH"]
    dat_full["DEPTH_ND"] = dat_full["DEPTH_AM"] - dat_full["DEPTH"]

    # If all ND depths are 0 → set VAF_ND = 0 and warn
    if (dat_full["DEPTH_ND"] == 0).all():
        print("ALT_DEPTH_AM is the same as ALT_DEPTH, you should not use VAF_ND")
        dat_full["VAF_ND"] = 0
    else:
        dat_full["VAF_ND"] = np.where(
            dat_full["DEPTH_ND"] > 0,
            dat_full["ALT_DEPTH_ND"] / dat_full["DEPTH_ND"],
            0
        )


    # compute VAF distortion
    dat_full["VAF_distortion"] = dat_full["VAF_AM"] / dat_full["VAF"]
    dat_full["VAF_distortion_sq"] = np.log10(dat_full["VAF"]) / np.log10(dat_full["VAF_AM"])

    dat_full["VAF_distorted_expanded"] = dat_full["VAF_distortion"] > vaf_distortion_threshold
    dat_full["VAF_distorted_expanded_sq"] = dat_full["VAF"] < ( dat_full["VAF_AM"] ** 1.5 )

    dat_full["VAF_distorted_expanded"] = dat_full["VAF_distorted_expanded"].fillna(True)
    dat_full["VAF_distorted_expanded_sq"] = dat_full["VAF_distorted_expanded_sq"].fillna(True)

    
    # define mean position in read 
    dat_full["PMEAN"] = dat_full["PMN"].astype(float) if "PMN" in dat_full.columns else -1
    dat_full["PSTD"] = dat_full["PST"].astype(float) if "PST" in dat_full.columns else -1


    # subset dataframe to the columns of interest
    dat_full_reduced = dat_full[ columns_to_keep ].reset_index(drop = True)


    # we generate an ID field that can be used for:
    #     - retrieving annotations
    #     - computing intersections between different mutation sets
    if list(dat_full_reduced.columns[:4]) == ['CHROM', 'POS', 'REF', 'ALT']:
        dat_full_reduced["MUT_ID"] = [ f'{row[0]}:{row[1]}_{row[2]}>{row[3]}' for _, row in dat_full_reduced.iterrows() ]
        print("MUT_ID generated")
    else:
        print("No ID field generated, 'CHROM', 'POS', 'REF', 'ALT' not in the first positions of the df")


    return dat_full_reduced

def reformat_alleles_to_OncodriveFML(dat, letters = ['A', 'T', 'C', 'G']):
    """
    This function receives a ROW of a dataframe
    with the mutation coordinates produced by VarDictJava
    at least these columns should be there:
    ['CHROM', 'POS', 'REF', 'ALT', 'Location', 'Allele']

    'CHROM', 'POS', 'REF', 'ALT'
    come from the original vardict file

    'Location', 'Allele'
    are added when merging the VEP annotation

    the function returns 4 values correponding to the updated coordinates, and alleles
    corresponding to the format required by tools such as OncodriveFML
    (indels represented by - without repeating the nucleotides.)
    """
    chromosome = dat["CHROM"].replace("chr", "")

    if dat["ALT"][0] in letters:
        position = dat["Location"].split(":")[1].split("-")[0]
        if len(dat["REF"]) == len(dat["ALT"]):
            allele = dat["Allele"]
            ref = dat["REF"]

        elif len(dat["REF"]) > len(dat["ALT"]):
            allele = dat["Allele"]
            ref = dat["REF"][len(dat["ALT"]):]
#            TTACGA     T
#            TACGA

        else:
            allele = dat["Allele"]
            ref = "-"

        return chromosome, position, ref, allele

    return chromosome, dat['POS'], dat['REF'], dat['ALT']

def add_alternative_format_columns(dataframe):
    """
    This function receives a dataframe
    with the mutation coordinates produced by VarDictJava
    at least these columns should be there:
    ['CHROM', 'POS', 'REF', 'ALT', 'Location', 'Allele']

    'CHROM', 'POS', 'REF', 'ALT'
    come from the original vardict file

    'Location', 'Allele'
    are added when merging the VEP annotation


    """
    for col in ['CHROM', 'POS', 'REF', 'ALT', 'Location', 'Allele']:
        c = dataframe[col]  # Trying to access a non-existent column

    # here we compute a new variant format that will be useful for running OncodriveFML
    new_format_muts = pd.DataFrame.from_records(
        dataframe[['CHROM', 'POS', 'REF', 'ALT', 'Location', 'Allele']].apply(
            reformat_alleles_to_OncodriveFML, axis = 1
        ).values
    )

    new_format_muts.columns = ["CHROM_ensembl", "POS_ensembl", "REF_ensembl", "ALT_ensembl"]

    return pd.concat((dataframe, new_format_muts), axis = 1)

def update_indel_info(df):
    """
    Update the DataFrame with INDEL_LENGTH and INDEL_INFRAME for indel rows,
    while setting default values for other rows.

    Parameters:
    df (pd.DataFrame): Original DataFrame containing the mutation data.

    Returns:
    pd.DataFrame: Updated DataFrame with INDEL_LENGTH and INDEL_INFRAME columns.
    """

    # Ensure the DataFrame has the necessary columns
    required_columns = ["TYPE", "REF", "ALT", "canonical_Consequence_broader"]
    for col in required_columns:
        if col not in df.columns:
            raise ValueError(f"DataFrame is missing required column: {col}")

    # Filter to get indices of rows that are indels (INSERTION or DELETION)
    indel_positions = df.loc[df["TYPE"].isin(["INSERTION", "DELETION"])].index

    # Initialize new columns with default values
    df['INDEL_LENGTH'] = 0
    df['INDEL_INFRAME'] = False
    df['INDEL_MULTIPLE3'] = False

    # Create a DataFrame specifically for indels
    indel_maf_df = df.loc[indel_positions].copy()

    # Calculate INDEL_LENGTH
    indel_maf_df["INDEL_LENGTH"] = (indel_maf_df["REF"].str.len() - indel_maf_df["ALT"].str.len()).abs()

    # Determine if the indel is in-frame
    indel_maf_df["INDEL_MULTIPLE3"] = indel_maf_df["INDEL_LENGTH"] % 3 == 0
    indel_maf_df["INDEL_INFRAME"] = indel_maf_df["INDEL_LENGTH"] % 3 == 0

    # Apply additional conditions to INDEL_INFRAME
    indel_maf_df.loc[indel_maf_df["INDEL_LENGTH"] >= 15, "INDEL_INFRAME"] = False
    indel_maf_df.loc[indel_maf_df["canonical_Consequence_broader"] == 'nonsense', "INDEL_INFRAME"] = False
    indel_maf_df.loc[indel_maf_df["canonical_Consequence_broader"] == 'essential_splice', "INDEL_INFRAME"] = False

    # Update the original DataFrame with the calculated values for indels
    df.loc[indel_positions, 'INDEL_LENGTH'] = indel_maf_df['INDEL_LENGTH']
    df.loc[indel_positions, 'INDEL_INFRAME'] = indel_maf_df['INDEL_INFRAME']

    return df

def vartype(x,
            letters = ['A', 'T', 'C', 'G'],
            len_SV_lim = 50
            ):
    """
    Define the TYPE of a variant
    """
    if ">" in (x["REF"] + x["ALT"]) or "<" in (x["REF"] + x["ALT"]):
        return "SV"

    elif len(x["REF"]) > (len_SV_lim+1) or len(x["ALT"]) > (len_SV_lim+1) :
        return "SV"

    elif x["REF"] in letters and x["ALT"] in letters:
        return "SNV"

    elif len(x["REF"]) == len(x["ALT"]):
        return "MNV"

    elif x["REF"] == "-" or ( len(x["REF"]) == 1 and x["ALT"].startswith(x["REF"]) ):
        return "INSERTION"

    elif x["ALT"] == "-" or ( len(x["ALT"]) == 1 and x["REF"].startswith(x["ALT"]) ):
        return "DELETION"

    return "COMPLEX"


def main_vcf_reader(vcf, sampleid, level, vaf_distortion_threshold):
    keep_all_columns = [
        "CHROM", "POS", "REF", "ALT", "FILTER", "INFO", "FORMAT",
        "SAMPLE", "DEPTH", "ALT_DEPTH", "REF_DEPTH", "VAF",
        'vd_DEPTH', 'vd_ALT_DEPTH', 'vd_REF_DEPTH', 'vd_VAF',
        "numNs", 'VAF_Ns',
        'DEPTH_AM', 'ALT_DEPTH_AM', 'REF_DEPTH_AM', "VAF_AM",
        "numNs_AM", "VAF_Ns_AM",
        'DEPTH_ND', 'ALT_DEPTH_ND',
        "VAF_ND",
        "VAF_distorted_expanded",
        "VAF_distorted_expanded_sq",
        "VAF_distortion",
        "VAF_distortion_sq",
        "PMEAN", "PSTD"
    ]

    sample_muts = read_from_vardict_VCF_all(
        sampleid,
        vcf,
        columns_to_keep=keep_all_columns,
        vaf_distortion_threshold=vaf_distortion_threshold
    )

    sample_muts["SAMPLE_ID"] = sampleid
    sample_muts["METHOD"] = level

    # remove non canonical chromosomes
    sample_muts = sample_muts[sample_muts["CHROM"].apply(lambda x: '_' not in x)].reset_index(drop=True)


    sample_muts["TYPE"] = sample_muts[["REF", "ALT"]].apply(vartype, axis=1)

    return sample_muts




if __name__ == '__main__':
    main()
