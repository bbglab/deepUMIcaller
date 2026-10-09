#!/usr/bin/env python
"""
Estimate UMI collision probabilities for a single sample.

Inputs (produced by the deepUMIcaller pipeline):
  - one or more `*.position_group_sizes.txt` files from FGUMI_GROUPREADSBYUMI
    (fgumi group step), describing how many fragments share each position group
  - the `*.umi_counts.txt` file from FGUMI_COLLECTDUPLEXSEQMETRICS
    (fgumi duplex-metrics step), with the distribution of single UMI observations
  - the `*.duplex_umi_counts.txt` file from FGUMI_COLLECTDUPLEXSEQMETRICS,
    with the distribution of duplex UMI observations

Output:
  - a TSV file with the UMI collision summary statistics for the sample
  - optionally, a PDF plot of the UMI frequency distributions
"""

import os
import sys

# Use a non-interactive backend: this runs inside a container without a display
import matplotlib
matplotlib.use("Agg")

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import click
from scipy.optimize import brentq


# --------------------------------------------------------------------------
# Input parsing
# --------------------------------------------------------------------------

def read_position_group_sizes(position_group_sizes_files):
    """
    Read one or more `*.position_group_sizes.txt` files (fgumi group step) and
    combine them into a single dataframe.

    position_group_sizes_files: list of paths
    """
    file_list = [position_group_sizes_files] if isinstance(position_group_sizes_files, (str, os.PathLike)) else list(position_group_sizes_files)
    sample_counts = pd.concat([pd.read_table(f) for f in file_list], ignore_index=True)
    sample_counts = sample_counts.groupby(["position_group_size"], as_index=False).agg({"count": "sum"})
    sample_counts["total_fragments"] = sample_counts["position_group_size"] * sample_counts["count"]
    sample_counts["repeated"] = sample_counts["position_group_size"] > 1
    return sample_counts


def get_single_umi_probabilities_vector(umi_counts_files, verbose=False):
    """
    Read the `*.umi_counts.txt` file (fgumi duplex-metrics step) and extract the
    probabilities of the main (high-frequency) UMI tags.

    Returns:
        unique_tags_number: number of main UMI tags (elbow point of the distribution)
        umi_tag_probabilities: raw fractions of the main UMI tags
        umi_tag_probabilities_norm: normalized probabilities of the main UMI tags
    """
    file_list = [umi_counts_files] if isinstance(umi_counts_files, (str, os.PathLike)) else list(umi_counts_files)
    umi_counts_data = pd.concat([pd.read_table(f) for f in file_list], ignore_index=True)
    umi_counts_data = umi_counts_data.groupby(["umi"], as_index=False).agg({"unique_observations": "sum"})
    umi_counts_data["fraction_unique_observations"] = umi_counts_data["unique_observations"] / umi_counts_data["unique_observations"].sum()
    umi_counts_data = umi_counts_data.sort_values("fraction_unique_observations", ascending=False).reset_index(drop=True)
    umi_counts_data["cummulative_fraction"] = umi_counts_data["fraction_unique_observations"].cumsum()

    # identify the elbow point in the cumulative distribution of UMIs to determine a cutoff for high-frequency UMIs
    init_val = 1
    i = 0
    diff = 0
    while diff < (init_val / 2):
        new_val = umi_counts_data.loc[i, "fraction_unique_observations"]
        if verbose:
            print(i, new_val, init_val)
        if i == 0:
            init_val = umi_counts_data.loc[i, "fraction_unique_observations"]
        else:
            diff = init_val - umi_counts_data.loc[i, "fraction_unique_observations"]
        init_val = umi_counts_data.loc[i, "fraction_unique_observations"]
        if verbose:
            print(diff)
        i += 1

    umi_tag_probabilities = umi_counts_data.loc[:i - 1, "fraction_unique_observations"].values
    umi_tag_probabilities_norm = umi_tag_probabilities / umi_tag_probabilities.sum()

    return i, umi_tag_probabilities, umi_tag_probabilities_norm


def get_duplex_umi_probabilities_vector(duplex_umi_counts_files, size, verbose=False):
    """
    Read the `*.duplex_umi_counts.txt` file (fgumi duplex-metrics step) and extract
    the probabilities of the most frequent duplex UMI tags.

    Returns:
        duplex_tag_probabilities: raw fractions of the top `size` duplex UMI tags
        duplex_tag_probabilities_norm: normalized probabilities
    """
    file_list = [duplex_umi_counts_files] if isinstance(duplex_umi_counts_files, (str, os.PathLike)) else list(duplex_umi_counts_files)
    umi_counts_data = pd.concat([pd.read_table(f) for f in file_list], ignore_index=True)
    print(umi_counts_data.head())
    umi_counts_data = umi_counts_data.groupby(["umi"], as_index=False).agg({"unique_observations": "sum"})
    umi_counts_data["fraction_unique_observations"] = umi_counts_data["unique_observations"] / umi_counts_data["unique_observations"].sum()
    umi_counts_data = umi_counts_data.sort_values("fraction_unique_observations", ascending=False).reset_index(drop=True)
    umi_counts_data["cummulative_fraction"] = umi_counts_data["fraction_unique_observations"].cumsum()

    umi_tag_probabilities = umi_counts_data.loc[:size, "fraction_unique_observations"].values
    umi_tag_probabilities_norm = umi_tag_probabilities / umi_tag_probabilities.sum()

    return umi_tag_probabilities, umi_tag_probabilities_norm


# --------------------------------------------------------------------------
# Functions for the estimation of UMI collisions
# --------------------------------------------------------------------------

def posterior_expectation_one_repeat(p):
    """
    p: numpy array of probabilities
    """
    num = np.sum(p / ((1 - p) ** 2))
    den = np.sum(p / (1 - p))
    return num / den


def func(N, p):
    """
    p: numpy array of probabilities
    N: number of independent trials
    """
    return np.sum(1 - (1 - p) ** N)


def find_N(p, target):
    """
    p: numpy array of probabilities
    target: target value for the sum
    """
    # Define a function that we want to find the root of
    def objective(N):
        return func(N, p) - target

    # Use brentq to find the root of the objective function
    N_solution = brentq(objective, 1, 1e6)  # Search for N in the range [1, 1e6]

    return N_solution


def compute_original_fragments_per_sample(counts_df, vector_of_probabilities):
    """
    For each observed position group size (number of fragments sharing a cut site),
    estimate the original number of independent fragments that gave rise to it.
    """
    n_repeats_observed = counts_df["position_group_size"].unique().tolist()
    if 1 in n_repeats_observed:
        n_repeats_observed.remove(1)

    mapping_observed_repeats_to_original_cuts = {1: posterior_expectation_one_repeat(vector_of_probabilities).item()}

    for n_repeats in n_repeats_observed:
        N_solution = find_N(vector_of_probabilities, n_repeats)
        mapping_observed_repeats_to_original_cuts[n_repeats] = N_solution

    return mapping_observed_repeats_to_original_cuts


def plot_umi_frequencies(x, type_of_umi="duplex", sample="sample", output_file=None):
    data = pd.DataFrame(sorted(x, reverse=True))
    data.columns = ["frequency"]
    data = data.reset_index()

    # plot cumulative distribution of UMIs (fraction_unique_observations)
    sns.lineplot(data=data, x="index", y="frequency")
    plt.xlabel(f"{type_of_umi.capitalize()} UMI Rank")
    plt.ylabel("Fraction of Raw Observations")
    plt.title(f"{sample}\nCumulative Distribution of {type_of_umi.capitalize()} UMIs")
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


# --------------------------------------------------------------------------
# Main computation for a single sample
# --------------------------------------------------------------------------

def compute_umi_collisions(sample_name, position_group_sizes_files, umi_counts_files,
                           duplex_umi_counts_files, output_file, plot_file=None, verbose=False):
    """
    Compute the UMI collision statistics for a single sample and write them to a TSV file.
    """
    unique_tags_number, umi_tag_probabilities, umi_tag_probabilities_norm = \
        get_single_umi_probabilities_vector(umi_counts_files, verbose=verbose)

    print(f"{sample_name}: main_UMI_tags\t= {unique_tags_number}")

    # assume the two sides of the duplex are drawn independently from the same UMI distribution
    tag_probabilities_assume_random_side = np.matmul(
        umi_tag_probabilities_norm.reshape(-1, 1),
        umi_tag_probabilities_norm.reshape(1, -1)
    ).flatten()

    duplex_tag_probabilities, duplex_tag_probabilities_norm = \
        get_duplex_umi_probabilities_vector(duplex_umi_counts_files, size=unique_tags_number ** 2, verbose=verbose)

    vector_of_probabilities = tag_probabilities_assume_random_side

    if plot_file:
        plot_umi_frequencies(umi_tag_probabilities_norm, type_of_umi="single", sample=sample_name,
                             output_file=plot_file)

    unique_tags = unique_tags_number
    prob_vector = vector_of_probabilities

    probability_of_same_tag_2_frag = np.sum(prob_vector ** 2)
    D_eff = 1 / probability_of_same_tag_2_frag
    print(f"{sample_name}: D_eff\t\t= {D_eff:.0f}")

    sample_count_df = read_position_group_sizes(position_group_sizes_files)

    try:
        mapping_observed_repeats_to_original_cuts = compute_original_fragments_per_sample(sample_count_df, prob_vector)
    except Exception as e:
        print(f"Error computing original fragments for {sample_name}: {e}")
        summary_stats_df = pd.DataFrame(columns=["sample", "unique_tags_number", "unique_tags_number_duplex", "D_eff", "lost_fragments", "lost_proportion"])
        summary_stats_df.to_csv(output_file, sep="\t", index=False)
        print(f"{sample_name}: summary written to {output_file}")
        return None

    sample_count_df["original_fragments_per_cutsite"] = sample_count_df["position_group_size"].map(mapping_observed_repeats_to_original_cuts)
    sample_count_df["original_fragments"] = sample_count_df["original_fragments_per_cutsite"] * sample_count_df["count"]
    lost_fragments = sample_count_df["original_fragments"].sum() - sample_count_df["total_fragments"].sum()

    print(f"{sample_name}: Lost fragments\t= {lost_fragments:,.0f}")

    lost_proportion = lost_fragments / sample_count_df["original_fragments"].sum()

    print(f"{sample_name}: Lost proportion\t= {lost_proportion:.2%}")

    summary_stats_df = pd.DataFrame(
        [[sample_name, unique_tags, unique_tags ** 2, round(D_eff), round(lost_fragments), lost_proportion]],
        columns=["sample", "unique_tags_number", "unique_tags_number_duplex", "D_eff",
                 "lost_fragments", "lost_proportion"]
    )
    summary_stats_df["percentage_lost"] = summary_stats_df["lost_proportion"] * 100
    summary_stats_df["lost_proportion"] = summary_stats_df["lost_proportion"].round(7)

    summary_stats_df.to_csv(output_file, sep="\t", index=False)
    print(f"{sample_name}: summary written to {output_file}")

    return summary_stats_df


@click.command()
@click.option('--sample-name', '-n', required=True, type=str, help='Name of the sample.')
@click.option('--position-group-sizes', '-g', required=True, multiple=True,
              type=click.Path(exists=True),
              help='Path to the position_group_sizes.txt file(s) from the fgumi group step. '
                   'Can be specified multiple times for aggregation (e.g. per-chromosome files).')
@click.option('--umi-counts', '-u', required=True, multiple=True,
              type=click.Path(exists=True),
              help='Path to the umi_counts.txt file from the fgumi duplex-metrics step.')
@click.option('--duplex-umi-counts', '-d', required=True, multiple=True,
              type=click.Path(exists=True),
              help='Path to the duplex_umi_counts.txt file from the fgumi duplex-metrics step.')
@click.option('--output-file', '-o', required=True, type=click.Path(),
              help='Path to the output TSV file with the UMI collision statistics.')
@click.option('--plot-file', '-p', type=click.Path(), default=None,
              help='Optional path to save a PDF plot of the UMI frequency distributions.')
@click.option('--verbose', '-v', is_flag=True, default=False, help='Print verbose output.')
def main(sample_name, position_group_sizes, umi_counts, duplex_umi_counts, output_file, plot_file, verbose):
    """
    Compute the UMI collision probabilities for a single sample from the fgumi group
    and fgumi duplex-metrics outputs.
    """
    compute_umi_collisions(
        sample_name=sample_name,
        position_group_sizes_files=list(position_group_sizes),
        umi_counts_files=list(umi_counts),
        duplex_umi_counts_files=list(duplex_umi_counts),
        output_file=output_file,
        plot_file=plot_file,
        verbose=verbose,
    )


if __name__ == '__main__':
    main()
