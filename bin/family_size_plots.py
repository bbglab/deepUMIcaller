#!/usr/bin/env python

import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


def compute_duplicates(file,
                        count_var = "count"):

    groupby_umi_families = pd.read_table(file, sep = "\t", header = 0)
    groupby_umi_families['reads'] = groupby_umi_families[count_var] * groupby_umi_families['family_size']

    unique_molecules = groupby_umi_families[count_var].sum()
    total_molecules = groupby_umi_families['reads'].sum()
    duplicated_molecules = total_molecules - unique_molecules

    return duplicated_molecules / total_molecules    



def wrap_family_size_curves(data_sscs):
    '''
    wrap_family_size_curves(data_scss_grouped[["family_size", "in_duplex", "fraction_reads"]])
    '''

    # Pivot the DataFrame
    pivot_df = data_sscs.pivot_table(index='in_duplex', columns='family_size', values='fraction_reads', fill_value=0)

    # If necessary, add columns for family sizes that are not present in the data
    missing_cols = []
    for i in range(1, 202):
        if i not in pivot_df.columns:
            missing_cols.append(i)
    pivot_df[missing_cols] = 0

    # Sort the columns
    pivot_df = pivot_df.reindex(sorted(pivot_df.columns), axis=1)

    return pivot_df.reset_index()


def compute_family_sizes_curve(sample, duplex_fam_data, prefix_figure, percent_duplicates, confidence = '2 1 1', confidence_name = 'low'):
    '''
    this function computes the family size metrics curve
    '''
    confidence_list = [int(x) for x in confidence.split(' ')]
    min_reads = confidence_list[0]
    strand1 = confidence_list[1]
    strand2 = confidence_list[2]

    duplex_fam_data["in_duplex"] = np.logical_and( (duplex_fam_data["ab_size"] + duplex_fam_data["ba_size"]) >= min_reads,
                                                        np.logical_and(duplex_fam_data["ab_size"] >= strand1,
                                                                        duplex_fam_data["ba_size"] >= strand2
                                                                        )
                                                    )

    data_duplex_families_small = duplex_fam_data[["ab_size", "ba_size", "count", "in_duplex"]].copy()

    family_size1 = data_duplex_families_small[["ab_size", "count", "in_duplex"]]
    family_size2 = data_duplex_families_small[["ba_size", "count", "in_duplex"]]
    family_size1.columns = ["family_size", "count", "in_duplex"]
    family_size2.columns = ["family_size", "count", "in_duplex"]

    data_scss = pd.concat((family_size1, family_size2))
    data_scss = data_scss[data_scss["family_size"] > 0].reset_index(drop = True)
    data_scss["count_reads"] = data_scss["family_size"] * data_scss["count"]

    # this way we force that everything above 200 copies gets compressed into 201
    data_scss['family_size'] = data_scss['family_size'].apply(lambda x: x if x <= 200 else 201)
    
    data_scss_grouped = data_scss.groupby(["family_size", "in_duplex"]).sum().reset_index()
    data_scss_grouped["family_size"] = data_scss_grouped["family_size"].astype(int)
    data_scss_grouped["fraction"] = data_scss_grouped["count"] / data_scss_grouped["count"].sum()
    data_scss_grouped["fraction_reads"] = data_scss_grouped["count_reads"] / data_scss_grouped["count_reads"].sum()

    data_family_size_curve = wrap_family_size_curves(data_scss_grouped[["family_size", "in_duplex", "fraction_reads"]])
    orig_cols = list(data_family_size_curve.columns)
    data_family_size_curve["sample"] = sample
    data_family_size_curve["quality"] = confidence_name
    data_family_size_curve = data_family_size_curve[["sample",  "quality"] + orig_cols]


    total_duplex = data_duplex_families_small["count"][data_duplex_families_small["in_duplex"]].sum()
    total_scss = data_scss["count"].sum()
    total_reads = data_scss_grouped["count_reads"].sum()

    total_scss_duplex = data_scss["count"][data_scss["in_duplex"]].sum()
    total_reads_duplex = data_scss_grouped["count_reads"][data_scss["in_duplex"]].sum()
    
    total_scss_nonduplex = data_scss["count"][~data_scss["in_duplex"]].sum()
    total_reads_nonduplex = data_scss_grouped["count_reads"][~data_scss["in_duplex"]].sum()

    expected_dcs = round(total_scss / 2)
    recovery_of_dcs = total_duplex / expected_dcs * 100
    unique_reads = total_duplex + total_scss_nonduplex

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (14, 6))
    sns.lineplot(data = data_scss_grouped,
                x = "family_size",
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
    ax2.axis('off')
    fig.suptitle(f"{sample} Duplex:{confidence_name}({confidence})")
    plt.show()
    fig.savefig(f"{prefix_figure}.{confidence_name}.pdf", bbox_inches='tight')

    try:
        max_indices = np.argmax( data_scss_grouped[data_scss_grouped['in_duplex'] == True]['fraction_reads'].values )
        peak_size = data_scss_grouped[data_scss_grouped['in_duplex'] == True]['family_size'].values[max_indices]
    except:
        peak_size = 0


    keys = ['sample', 'quality', 'raw_reads',
                'duplicates', 'sscs', 'dcs',
                'expected_dcs', 'recovery_of_dcs', 'unique_reads', 'uq_reads_duplex', 'unique_molecules',
                'raw_x_dcs', 'raw_x_sscs', 'sscs_x_dcs',
                'duplex_raw_reads', 'duplex_sscs', 'duplex_raw_x_dcs', 'duplex_raw_x_sscs',
                'noduplex_raw_reads', 'noduplex_sscs', 'noduplex_raw_x_sscs',
                'peak_size']

    values = [sample, confidence_name, total_reads,
                round(percent_duplicates, 3), total_scss,
                total_duplex,
                expected_dcs, round(recovery_of_dcs, 3), unique_reads, round(total_duplex/unique_reads*100, 3), round(unique_reads/2),
                round(total_reads/total_duplex, 3), round(total_reads/total_scss, 3), round(total_scss/total_duplex, 3),
                total_reads_duplex, total_scss_duplex, round(total_reads_duplex/total_duplex, 3), round(total_reads_duplex/total_scss_duplex, 3),
                total_reads_nonduplex, total_scss_nonduplex, round(total_reads_nonduplex/total_scss_nonduplex, 3),
                peak_size]

    sample_data = pd.DataFrame([values], columns = keys)

    return sample_data, data_family_size_curve



def stats_fam_size2plot(sample, groupby_metrics_file, duplex_metrics_file, output_prefix):
    
    # compute duplicate rate from groupby stats
    prop_duplicates = compute_duplicates(groupby_metrics_file)
    percent_duplicates = prop_duplicates * 100

    # compute family size distributions from duplex stats data
    data_duplex_families = pd.read_table(f"{duplex_metrics_file}")


    sample_data_high, data_family_size_curve_high = compute_family_sizes_curve(sample, data_duplex_families,
                                                                                prefix_figure = output_prefix,
                                                                                percent_duplicates=percent_duplicates,
                                                                                confidence = '6 3 3', confidence_name = 'high')
    sample_data_med, data_family_size_curve_med = compute_family_sizes_curve(sample, data_duplex_families,
                                                                                    prefix_figure = output_prefix,
                                                                                    percent_duplicates=percent_duplicates,
                                                                                confidence = '4 2 2', confidence_name = 'med')
    sample_data_low, data_family_size_curve_low = compute_family_sizes_curve(sample, data_duplex_families, 
                                                                                    prefix_figure = output_prefix,
                                                                                    percent_duplicates=percent_duplicates,
                                                                                    confidence = '2 1 1', confidence_name = 'low')

    sample_data = pd.concat((sample_data_high, sample_data_med, sample_data_low), axis = 0)
    sample_data.to_csv(f"{output_prefix}.sample_data.tsv",
                                    sep = "\t",
                                    header = True,
                                    index = False)


    family_size_curve = pd.concat((data_family_size_curve_high, data_family_size_curve_med, data_family_size_curve_low), axis = 0)
    family_size_curve.to_csv(f"{output_prefix}.family_curve.tsv",
                                    sep = "\t",
                                    header = True,
                                    index = False)

    return 0


sam = sys.argv[1]
groupby_metrics_file = sys.argv[2]
duplex_metrics_file = sys.argv[3]
output_file = sys.argv[4]

global limx
limx = (0,50)

stats_fam_size2plot(sam, groupby_metrics_file, duplex_metrics_file, output_file)