#!/usr/bin/env python

import click
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

# TODO
# pending to add some parsing of the confidence level to the outputs

mutation_types = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
available_changes = mutation_types + ['all']

snv_color = {
        'C>A' : '#03bcee',
        'C>G' : '#010101',
        'C>T' : '#e32926',
        'T>A' : '#cac9c9',
        'T>C' : '#a1ce63',
        'T>G' : '#ebc6c4',
        'all' : '#FF00FF'
    }


def compute_ratios(df_valid, number_of_initial_positions):
    df_valid["all"] = df_valid[mutation_types].sum(axis = 1)
    ratios = {}
    for col in available_changes:
        values = df_valid[col].values
        first_positions = values[:number_of_initial_positions]
        rest = values[number_of_initial_positions:]
        mean_first_positions = first_positions.mean()
        mean_rest = rest.mean() if len(rest) > 0 else 0
        if mean_rest == 0:
            ratio = -1
        else:
            ratio = mean_first_positions / mean_rest
        ratios[col] = ratio

    return ratios



def plot_samples_lines(df, number_of_initial_positions, pdf = None):
    plt.figure(figsize=(12, 6))
    for sample in df.index:
        plt.plot(df.columns, df.loc[sample], marker='o', label=sample)
    plt.axhline(2, color='red', linestyle='--', label='Ratio = 2')
    plt.axhline(1, color='grey', linestyle='--', label='Ideal scenario')
    plt.title(f'Mutation ratios (mean first {number_of_initial_positions} / mean rest) by mutation type, per sample')
    plt.ylabel('Ratio')
    plt.xlabel('Mutation type')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', frameon=False)
    plt.tight_layout()
    if pdf is not None:
        try:
            pdf.savefig()
        except Exception:
            plt.savefig('ratios_samples_lines.png')
    else:
        plt.savefig('ratios_samples_lines.png')
    plt.close()

def plot_ratios_lines(df, number_of_initial_positions, pdf = None):
    plt.figure(figsize=(12, 6))
    for col in df.columns:
        color = snv_color.get(col, 'grey')
        plt.plot(df.index, df[col], marker='o', label=col, color=color)
    plt.axhline(2, color='red', linestyle='--', label='Ratio = 2')
    plt.axhline(1, color='grey', linestyle='--', label='Ideal scenario')
    plt.title(f'Mutation ratios (mean first {number_of_initial_positions} / mean rest) by sample')
    plt.ylabel('Ratio')
    plt.xlabel('Sample')
    plt.xticks(rotation=90)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small', frameon=False)
    plt.tight_layout()
    if pdf is not None:
        try:
            pdf.savefig()
        except Exception:
            plt.savefig('ratios_lines.png')
    else:
        plt.savefig('ratios_lines.png')
    plt.close()


# Plotting function for heatmaps
def plot_ratios_heatmap(df, number_of_initial_positions, pdf = None):
    plt.figure(figsize=(12, max(6, 0.4*len(df))))
    sns.heatmap(df, annot=True, fmt=".2f", cmap="coolwarm", center=1, cbar_kws={'label': 'Ratio'})
    plt.title(f'Heatmap of mutation ratios (mean first {number_of_initial_positions} / mean rest)')
    plt.ylabel('Sample')
    plt.xlabel('Mutation type')
    plt.tight_layout()
    if pdf is not None:
        try:
            pdf.savefig()
        except Exception:
            plt.savefig('ratios_heatmap.png')
    else:
        plt.savefig('ratios_heatmap.png')
    plt.close()

def plot_logfc_heatmap(df, number_of_initial_positions, pdf = None):
    logfc_df = np.log2(df)
    plt.figure(figsize=(12, max(6, 0.4*len(df))))
    sns.heatmap(logfc_df, annot=True, fmt=".2f", cmap="coolwarm", vmin=-2, vmax=2, cbar_kws={'label': 'log2(fold change)', 'ticks': [-2, -1, 0, 1, 2]})
    plt.title(f'Heatmap of log2 fold change (mean first {number_of_initial_positions} / mean rest)')
    plt.ylabel('Sample')
    plt.xlabel('Mutation type')
    plt.tight_layout()
    if pdf is not None:
        try:
            pdf.savefig()
        except Exception:
            plt.savefig('ratios_logfc_heatmap.png')
    else:
        plt.savefig('ratios_logfc_heatmap.png')
    plt.close()



def get_sample_name(file_path):
    # Extract the filename
    filename = os.path.basename(file_path)
    # Find the part before 'MutsPer'
    match = re.match(r'(.+?)MutsPer', filename)
    if match:
        return match.group(1).rstrip('_-.')
    else:
        return filename




# Use click for command-line interface
@click.command()
@click.option('--initial', default=5, show_default=True, help='Number of initial positions to use')
def main(initial):
    number_of_initial_positions = initial

    # Recursively find all *MutsPerCycle.dat.csv files
    csv_files = glob.glob(os.path.join('./', '**', '*MutsPerCycle_mqc.csv'), recursive=True)

    ratios_per_file = {}

    for file_path in csv_files:
        sample_name = get_sample_name(file_path)
        print(f"Processing file: {file_path} (Sample: {sample_name})")
        try:
            df = pd.read_csv(file_path)
            df_valid = df[df.loc[:, "Count"] != 0]
            df_valid = df_valid.iloc[4:, :]
            if len(df_valid) <= 15:
                print("not enough values")
                continue
            ratios = compute_ratios(df_valid, number_of_initial_positions)

            ratios_per_file[sample_name] = ratios
        except Exception as e:
            print(f"Error processing {file_path}: {e}")

    ratios_df = pd.DataFrame.from_dict(ratios_per_file, orient='index')
    ratios_df.index.name = 'Sample'
    if 'all' in ratios_df.columns:
        ratios_df = ratios_df.sort_values(by='all', ascending=False)
    ratios_df.to_csv('ratios_per_sample_mqc.tsv', sep='\t')
    print("\nSaved ratios table to ratios_per_sample_mqc.tsv (sorted by 'all' ratio descending)")

    with PdfPages('mutation_ratios_summary_mqc.pdf') as pdf:
        plot_ratios_heatmap(ratios_df, number_of_initial_positions, pdf)
        plot_logfc_heatmap(ratios_df, number_of_initial_positions, pdf)
        plot_ratios_lines(ratios_df, number_of_initial_positions, pdf)
        plot_samples_lines(ratios_df, number_of_initial_positions, pdf)
    print("Saved all plots to mutation_ratios_summary_mqc.pdf")

    print(f"Ratio of mean(first {number_of_initial_positions} positions) / mean(rest):")
    print(ratios_df)

    print("\nSamples where the ratio > 2 for any column:")
    passing_samples = []
    for sample, ratios in ratios_per_file.items():
        passing_types = [mut_type for mut_type, ratio in ratios.items() if ratio > 2]
        if passing_types:
            print(f"  {sample} in ", end="")
            print(", ".join(passing_types))
            passing_samples.append({"Sample": sample, "MutationTypes": ", ".join(passing_types)})

    # Write summary of samples passing the threshold to a file
    if passing_samples:
        passing_df = pd.DataFrame(passing_samples)
        passing_df.to_csv('samples_passing_ratio_threshold_mqc.tsv', sep='\t', index=False)
        print("\nSaved summary of samples passing the threshold to samples_passing_ratio_threshold_mqc.tsv")
    else:
        print("\nNo samples passed the ratio > 2 threshold.")

if __name__ == '__main__':
    main()

