#!/usr/bin/env python

import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def compute_duplicates(file,
                        count_var = "count"):

    groupby_umi_families = pd.read_table(file, sep = "\\t", header = 0)
    groupby_umi_families['reads'] = groupby_umi_families[count_var] * groupby_umi_families['family_size']

    unique_molecules = groupby_umi_families[count_var].sum()
    total_molecules = groupby_umi_families['reads'].sum()
    duplicated_molecules = total_molecules - unique_molecules

    return duplicated_molecules / total_molecules    

def stats_fam_size2plot(sample, groupby_metrics_file, duplex_metrics_file, limx):
    
    # compute duplicate rate from groupby stats
    prop_duplicates = compute_duplicates(groupby_metrics_file)
    percent_duplicates = prop_duplicates * 100

    # compute family size distributions from duplex stats data
    data_duplex_families = pd.read_table(f"{duplex_metrics_file}")

    data_duplex_families["in_duplex"] = np.logical_and(data_duplex_families["ab_size"] > 0,
                                                        data_duplex_families["ba_size"] > 0)

    data_duplex_families["in_duplex_high"] = np.logical_and(data_duplex_families["ab_size"] >= 3,
                                                            data_duplex_families["ba_size"] >= 3)

    data_duplex_families_small = data_duplex_families[["ab_size", "ba_size", "count", "in_duplex", "in_duplex_high"]]
    
    family_size1 = data_duplex_families_small[["ab_size", "count", "in_duplex"]]
    family_size2 = data_duplex_families_small[["ba_size", "count", "in_duplex"]]
    family_size1.columns = ["family_size", "count", "in_duplex"]
    family_size2.columns = ["family_size", "count", "in_duplex"]

    data_scss = pd.concat((family_size1, family_size2))
    data_scss = data_scss[data_scss["family_size"] > 0].reset_index(drop = True)
    data_scss["count_reads"] = data_scss["family_size"] * data_scss["count"]

    data_scss_grouped = data_scss.groupby(["family_size", "in_duplex"]).sum().reset_index()
    data_scss_grouped["family_size"] = data_scss_grouped["family_size"].astype(int)
    data_scss_grouped["fraction"] = data_scss_grouped["count"] / data_scss_grouped["count"].sum()
    data_scss_grouped["fraction_reads"] = data_scss_grouped["count_reads"] / data_scss_grouped["count_reads"].sum()

    total_duplex = data_duplex_families_small["count"][data_duplex_families_small["in_duplex"]].sum()
    total_scss = data_scss["count"].sum()
    total_reads = data_scss_grouped["count_reads"].sum()

    total_scss_duplex = data_scss["count"][data_scss["in_duplex"]].sum()
    total_reads_duplex = data_scss_grouped["count_reads"][data_scss["in_duplex"]].sum()
    
    total_scss_nonduplex = data_scss["count"][~data_scss["in_duplex"]].sum()
    total_reads_nonduplex = data_scss_grouped["count_reads"][~data_scss["in_duplex"]].sum()
    
    total_non_duplex = data_duplex_families_small["count"][~data_duplex_families_small["in_duplex"]].sum()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (14, 6))
    sns.lineplot(data = data_scss_grouped,
                x = "family_size",
                # y = "count",
                y = "fraction_reads",
                hue = "in_duplex",
                marker = 'o',
                palette= {True : "k", False : "r"},
                ax = ax1
                )

    ax1.set_xlim(limx)
    ax1.set_xlabel("Family size")
    ax1.set_ylabel("Fraction of reads")
    ax1.legend(title = "in duplex\nfamilies")
    
    ax2.text(0.15, 0.8, f"Raw reads:   {total_reads:,}\nDuplicates:      {percent_duplicates:.1f}%\n\nSSCS:            {total_scss:,}" )
    ax2.text(0.15, 0.6, f"Raw/DCS:          {total_reads/total_duplex:.3f}\nRaw/SSCS:         {total_reads/total_scss:.3f}\nSSCS/Duplex:    {total_scss/total_duplex:.3f}")
    ax2.text(0.15, 0.25, f"IN DUPLEX\nRaw:                {total_reads_duplex:,} ({total_reads_duplex/(total_reads)*100:.1f}%)\nSSCs:               {total_scss_duplex:,}\nDuplex:            {total_duplex:,}\nRaw/DCS:         {total_reads_duplex/total_duplex:.3f}\nRaw/SSCS:        {total_reads_duplex/total_scss_duplex:.3f}\n")
    ax2.text(0.15, 0.05, f"NO DUPLEX\nRaw:                {total_reads_nonduplex:,} ({total_reads_nonduplex/(total_reads)*100:.1f}%)\nSSCs:               {total_scss_nonduplex:,} ({total_scss_nonduplex/(total_scss)*100:.1f}%)\nRaw/SSCS:        {total_reads_nonduplex/total_scss_nonduplex:.3f}")
    #ax2.text(0.1, 0.0, f"NO DUPLEX\nRaw:           {total_reads_nonduplex:,}\nSSCs:  Raw/SSCS:      {total_reads_nonduplex/total_scss_nonduplex:.3f}")
    ax2.axis('off')
    
    fig.suptitle(sample)    
    plt.show()

    duplex_peak = data_scss_grouped[data_scss_grouped['in_duplex'] == True]['']


    simple_counts_header = "raw_reads\tduplicates\tsscs\traw_x_dcs\traw_x_sscs\tsscs_x_dcs"
    in_duplex_header = "duplex_raw_reads\tduplex_sscs\tduplex_dcs\tduplex_raw_x_dcs\tduplex_raw_x_sscs"
    no_duplex_header = "noduplex_raw_reads\tnoduplex_sscs\tnoduplex_raw_x_sscs"
    header = f"{simple_counts_header}\t{in_duplex_header}\t{no_duplex_header}\n"
    
    simple_counts_features = f"{total_reads:}\t{percent_duplicates:.3}\t{total_scss:}\t{total_reads/total_duplex:.3f}\t{total_reads/total_scss:.3f}\t{total_scss/total_duplex:.3f}"
    in_duplex_features = f"{total_reads_duplex:}\t{total_scss_duplex:}\t{total_duplex:}\t{total_reads_duplex/total_duplex:.3f}\t{total_reads_duplex/total_scss_duplex:.3f}"
    no_duplex_features = f"{total_reads_nonduplex:}\t{total_scss_nonduplex:}\t{total_reads_nonduplex/total_scss_nonduplex:.3f}"
    features = f"{simple_counts_features}\t{in_duplex_features}\t{no_duplex_features}"

    print(header+features)
    
    return fig


sam = sys.argv[1]
groupby_metrics_file = sys.argv[2]
duplex_metrics_file = sys.argv[3]
output_file = sys.argv[4]

x_axis_limits = (0,50)

figure = stats_fam_size2plot(sam,
                            groupby_metrics_file,
                            duplex_metrics_file,
                            x_axis_limits)

figure.savefig(output_file, bbox_inches='tight')