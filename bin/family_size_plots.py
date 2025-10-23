#!/usr/bin/env python

import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import click


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


def compute_family_sizes_curve(sample, duplex_fam_data, prefix_figure, confidence = '2 1 1', confidence_name = 'low'):
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

    # Handle case where there are no reads (avoid division by zero)
    if total_reads == 0:
        percent_duplicates = 0
        expected_dscs = 0
        recovery_of_dscs = 0
    else:
        percent_duplicates = (1 - total_scss/total_reads) * 100
        expected_dscs = round(total_scss / 2)
        recovery_of_dscs = total_duplex / expected_dscs * 100 if expected_dscs > 0 else 0
    
    unique_reads = total_duplex + total_scss_nonduplex

    try:
        max_indices = np.argmax( data_scss_grouped[data_scss_grouped['in_duplex']]['fraction_reads'].values )
        peak_size = data_scss_grouped[data_scss_grouped['in_duplex']]['family_size'].values[max_indices]
    except:
        peak_size = 0

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    sns.lineplot(data=data_scss_grouped,
                 x="family_size",
                 y="fraction_reads",
                 hue="in_duplex",
                 marker='o',
                 palette={True: "k", False: "r"},
                 ax=ax1)

    # Add a vertical bar indicating the peak size
    ax1.axvline(x=peak_size, color="blue", linestyle="--", label=f"Peak:{peak_size}")

    ax1.set_xlim(limx)
    ax1.set_xlabel("Family size")
    ax1.set_ylabel("Fraction of reads")
    ax1.legend(title="in duplex\nfamilies")

    # Handle division by zero in text display
    raw_dscs_ratio = total_reads/total_duplex if total_duplex > 0 else 0
    raw_sscs_ratio = total_reads/total_scss if total_scss > 0 else 0
    sscs_dscs_ratio = total_scss/total_duplex if total_duplex > 0 else 0
    
    duplex_pct = total_reads_duplex/(total_reads)*100 if total_reads > 0 else 0
    duplex_raw_dscs_ratio = total_reads_duplex/total_duplex if total_duplex > 0 else 0
    duplex_raw_sscs_ratio = total_reads_duplex/total_scss_duplex if total_scss_duplex > 0 else 0
    
    nonduplex_pct = total_reads_nonduplex/(total_reads)*100 if total_reads > 0 else 0
    nonduplex_sscs_pct = total_scss_nonduplex/(total_scss)*100 if total_scss > 0 else 0
    nonduplex_raw_sscs_ratio = total_reads_nonduplex/total_scss_nonduplex if total_scss_nonduplex > 0 else 0

    ax2.text(0.15, 0.8, f"Raw reads:   {total_reads:,}\nDuplicates:      {percent_duplicates:.1f}%\n\nSSCs:            {total_scss:,}")
    ax2.text(0.15, 0.6, f"Raw/DSCs:          {raw_dscs_ratio:.3f}\nRaw/SSCs:         {raw_sscs_ratio:.3f}\nSSCs/DSCs:    {sscs_dscs_ratio:.3f}")
    ax2.text(0.15, 0.25, f"IN DUPLEX\nRaw:                {total_reads_duplex:,} ({duplex_pct:.1f}%)\nSSCs:               {total_scss_duplex:,}\nDSCs:            {total_duplex:,}\nRaw/DSCs:         {duplex_raw_dscs_ratio:.3f}\nRaw/SSCs:        {duplex_raw_sscs_ratio:.3f}\n")
    ax2.text(0.15, 0.05, f"NO DUPLEX\nRaw:                {total_reads_nonduplex:,} ({nonduplex_pct:.1f}%)\nSSCs:               {total_scss_nonduplex:,} ({nonduplex_sscs_pct:.1f}%)\nRaw/SSCs:        {nonduplex_raw_sscs_ratio:.3f}")
    ax2.axis('off')
    fig.suptitle(f"{sample} Duplex:{confidence_name}({confidence})")
    plt.show()
    fig.savefig(f"{prefix_figure}.{confidence_name}.pdf", bbox_inches='tight')



    keys = ['sample', 'quality', 'quality_threshold', 'raw_reads',
                'duplicates', 'sscs', 'dscs',
                'expected_dscs', 'recovery_of_dscs', 'unique_reads', 'uq_reads_duplex', 'unique_molecules',
                'raw_x_dscs', 'raw_x_sscs', 'sscs_x_dscs',
                'duplex_raw_reads', 'duplex_sscs', 'duplex_raw_x_dscs', 'duplex_raw_x_sscs',
                'noduplex_raw_reads', 'noduplex_sscs', 'noduplex_raw_x_sscs',
                'peak_size']

    # Calculate ratios with zero checks
    raw_x_dscs = round(total_reads/total_duplex, 3) if total_duplex > 0 else 0
    raw_x_sscs = round(total_reads/total_scss, 3) if total_scss > 0 else 0
    sscs_x_dscs = round(total_scss/total_duplex, 3) if total_duplex > 0 else 0
    uq_reads_duplex = round(total_duplex/unique_reads*100, 3) if unique_reads > 0 else 0
    unique_molecules = round(unique_reads/2)
    duplex_raw_x_dscs = round(total_reads_duplex/total_duplex, 3) if total_duplex > 0 else 0
    duplex_raw_x_sscs = round(total_reads_duplex/total_scss_duplex, 3) if total_scss_duplex > 0 else 0
    noduplex_raw_x_sscs = round(total_reads_nonduplex/total_scss_nonduplex, 3) if total_scss_nonduplex > 0 else 0

    values = [sample, confidence_name, confidence, total_reads,
                round(percent_duplicates, 3), total_scss,
                total_duplex,
                expected_dscs, round(recovery_of_dscs, 3), unique_reads, uq_reads_duplex, unique_molecules,
                raw_x_dscs, raw_x_sscs, sscs_x_dscs,
                total_reads_duplex, total_scss_duplex, duplex_raw_x_dscs, duplex_raw_x_sscs,
                total_reads_nonduplex, total_scss_nonduplex, noduplex_raw_x_sscs,
                peak_size]

    sample_data = pd.DataFrame([values], columns = keys)

    return sample_data, data_family_size_curve



def stats_fam_size2plot(sample, duplex_metrics_files, output_prefix, confidence = '4 2 2'):
    
    # compute family size distributions from duplex stats data
    # Handle single file or multiple files (list)
    if isinstance(duplex_metrics_files, list):
        # Multiple files: aggregate them by summing counts for each (ab_size, ba_size) combination
        dfs = []
        for metrics_file in duplex_metrics_files:
            try:
                df = pd.read_table(metrics_file)
                dfs.append(df)
            except Exception as e:
                print(f"Warning: Could not read {metrics_file}: {e}", file=sys.stderr)
                continue
        
        if not dfs:
            raise ValueError("No valid metrics files could be read")
        
        # Concatenate all dataframes
        data_duplex_families = pd.concat(dfs, ignore_index=True)
        
        # Sum counts for identical ab_size/ba_size combinations
        # Keep only ab_size, ba_size, and count columns for aggregation
        data_duplex_families = data_duplex_families[['ab_size', 'ba_size', 'count']]
        data_duplex_families = data_duplex_families.groupby(['ab_size', 'ba_size'], as_index=False)['count'].sum()
        
        # Note: fraction and fraction_gt_or_eq_size will be recalculated by the downstream analysis
        # so we don't need to preserve them during aggregation
    else:
        # Single file - read and keep only necessary columns
        df = pd.read_table(duplex_metrics_files)
        data_duplex_families = df[['ab_size', 'ba_size', 'count']].copy()


    sample_data_duplex, data_family_size_curve_duplex = compute_family_sizes_curve(sample, data_duplex_families,
                                                                                    prefix_figure = output_prefix,
                                                                                    confidence = confidence, confidence_name = 'duplex')
    sample_data_low, data_family_size_curve_low = compute_family_sizes_curve(sample, data_duplex_families, 
                                                                                    prefix_figure = output_prefix,
                                                                                    confidence = '2 1 1', confidence_name = 'allm')

    sample_data = pd.concat((sample_data_duplex, sample_data_low), axis = 0)
    sample_data.to_csv(f"{output_prefix}.sample_data.tsv",
                                    sep = "\t",
                                    header = True,
                                    index = False)


    family_size_curve = pd.concat((data_family_size_curve_duplex, data_family_size_curve_low), axis = 0)
    family_size_curve.to_csv(f"{output_prefix}.family_curve.tsv",
                                    sep = "\t",
                                    header = True,
                                    index = False)

    return 0


@click.command()
@click.option('--sample-name', '-n', required=True, type=str, help='Name of the sample.')
@click.option('--input-file', '-i', required=True, multiple=True, type=click.Path(exists=True), help='Path to the input file(s). Can be specified multiple times for aggregation.')
@click.option('--output-file', '-o', required=True, type=click.Path(), help='Path to the output file.')
@click.option('--confidence-level', '-c', default='4 2 2', show_default=True, type=str, help='Confidence level for the analysis.')
def main(sample_name, input_file, output_file, confidence_level):
    """
    Main function to process the input file(s) and generate plots.
    Accepts one or more input files for aggregation.
    """

    global limx
    limx = (0,50)

    # Convert tuple to list, or handle single file
    input_files = list(input_file) if len(input_file) > 1 else input_file[0]
    
    stats_fam_size2plot(sample_name, input_files, output_file, confidence=confidence_level)


if __name__ == '__main__':
    main()
