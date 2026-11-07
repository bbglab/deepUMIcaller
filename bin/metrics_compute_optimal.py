#!/usr/bin/env python

import sys
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import click


# Family size metrics

def fam_metrics2df(file, prefix = 'all'):
    
    data = pd.read_csv(file, sep = '\t')
    
    # Convert relevant columns to numeric, using errors='coerce' to handle non-numeric values
    numeric_columns = ['raw_reads', 'duplicates', 'sscs', 'dscs',
                        'expected_dscs', 'recovery_of_dscs', 'unique_reads', 'uq_reads_duplex',
                        'unique_molecules', 'raw_x_dscs', 'raw_x_sscs', 'sscs_x_dscs',
                        'duplex_raw_reads', 'duplex_sscs', 'duplex_raw_x_dscs',
                        'duplex_raw_x_sscs', 'noduplex_raw_reads', 'noduplex_sscs',
                        'noduplex_raw_x_sscs', 'peak_size']
    
    # Convert the specified columns to numeric
    for col in numeric_columns:
        if col in data.columns:  # Check if column exists before conversion
            data[col] = pd.to_numeric(data[col], errors='coerce')
    
    # transform multiple rows per sample to single row per sample
    pivot_df = data.pivot_table(index=['sample', 'raw_reads', 'sscs', 'raw_x_sscs', 'duplicates'], columns='quality', aggfunc='first')
    pivot_df.columns = [f'{col[1]}.{col[0]}' for col in pivot_df.columns]
    for prune in ["unique_molecules", "unique_reads"]:
        pivot_df[prune] = pivot_df[f"allm.{prune}"].copy()
        pivot_df = pivot_df.drop([f"allm.{prune}", f"duplex.{prune}"], axis = 'columns')
    data = pivot_df.reset_index()
    data.columns = ["sample"] + [f'FamMetrics>>{prefix}.{x}' for x in data.columns[1:] ]
    data["sample"] = data["sample"].astype(str)
    
    return data



def fam_curves2df(file, prefix = 'all'):
    """
    This function reads a TSV file containing family size distribution data curves and reshapes it into a wide format.
    TODO: store only the low (any duplex read) and the duplex confidence ones,
        the duplex confidence one should come from the confidence level defined
        by the user when running the pipeline and it should be recorded.
    """
    data = pd.read_csv(file, sep = '\t')

    summarized_data = data[['sample', 'quality', 'in_duplex']].copy()
    
    vals = data[[ x for x in list(data.columns) if x not in ['sample', 'quality', 'in_duplex'] ]]
    summarized_data["vector"] = vals.apply(lambda x: x.values, axis = 1)
    
    summarized_data_wide = summarized_data.pivot_table(index=['sample'], columns=['quality', 'in_duplex'], aggfunc='first')
    summarized_data_wide.columns = [f'FamCurves>>{prefix}.{x}.{y}.sscs' for _, x, y in summarized_data_wide.columns]
    summarized_data_wide = summarized_data_wide.reset_index()
    summarized_data_wide[f"FamCurves>>{prefix}.dist.sscs"] = summarized_data_wide[[f"FamCurves>>{prefix}.allm.False.sscs", f"FamCurves>>{prefix}.allm.True.sscs"]].sum(axis = 1)
    summarized_data_wide["sample"] = summarized_data_wide["sample"].astype(str)
    
    return summarized_data_wide

# Family size metrics - off target metrics

def compute_off_targets_metrics(data_family):
    """
    With the information from on target and total reads compute the familymetrics for off targets.
    """
    
    off_target_df = data_family[["sample"]].copy()
    
    for column in ['sscs', 'raw_reads', 'unique_molecules', 'unique_reads']:
        off_target_df[f'FamMetrics>>off_target.{column}'] = data_family[f'FamMetrics>>all.{column}'] - data_family[f'FamMetrics>>on_target.{column}']

    off_target_df['FamMetrics>>off_target.raw_x_sscs'] = off_target_df['FamMetrics>>off_target.raw_reads'] / off_target_df['FamMetrics>>off_target.sscs']

    column = 'dscs'

    for level in ["allm", "duplex"]:
        off_target_df[f'FamMetrics>>off_target.{level}.{column}'] = data_family[f'FamMetrics>>all.{level}.{column}'] - data_family[f'FamMetrics>>on_target.{level}.{column}']

        off_target_df[f'FamMetrics>>off_target.{level}.raw_x_dscs'] = off_target_df[f'FamMetrics>>off_target.raw_reads'] / off_target_df[f'FamMetrics>>off_target.{level}.dscs']
        off_target_df[f'FamMetrics>>off_target.{level}.sscs_x_dscs'] = off_target_df[f'FamMetrics>>off_target.sscs'] / off_target_df[f'FamMetrics>>off_target.{level}.dscs']

    for column in ['dist.sscs',
                'allm.False.sscs',
                'allm.True.sscs',
                'duplex.False.sscs',
                'duplex.True.sscs']:
        
        off_target_df[f'FamCurves>>off_target.{column}'] = data_family[f'FamCurves>>all.{column}'] * data_family[f'FamMetrics>>all.sscs'] \
                                                            - data_family[f'FamCurves>>on_target.{column}'] * data_family[f'FamMetrics>>on_target.sscs']

        off_target_df[f'FamCurves>>off_target.{column}'] = off_target_df[f'FamCurves>>off_target.{column}'] / off_target_df[f'FamMetrics>>off_target.sscs']
        
    return off_target_df    



# Compute FOR OPTIMAL metrics
def compute_desired_sequencing_metrics(family_data,
                                    desired_amplification = 18):
    """
    Here we are computing the metrics that are used for estimating the optimal sequencing.
    """
    
    family_data['FamMetrics>>on_target.unique_molecules_perc'] = family_data['FamMetrics>>on_target.unique_molecules'] / family_data['FamMetrics>>all.unique_molecules']
    family_data['FamMetrics>>off_target.unique_molecules_perc'] = family_data['FamMetrics>>off_target.unique_molecules'] / family_data['FamMetrics>>all.unique_molecules']
    
    family_data['FamMetrics>>amplification_bias.on/off'] = ( family_data['FamMetrics>>on_target.raw_reads'] / family_data['FamMetrics>>on_target.unique_molecules'] ) \
                                                                / ( family_data['FamMetrics>>off_target.raw_reads'] / family_data['FamMetrics>>off_target.unique_molecules'] )


    #UM * %on * desired_amplification = reads for optimal on target
    #UM * %off * (desired_amplification / amplification_bias.on/off )  = reads for optimal off target

    family_data['FamMetrics>>molecules_for_optimal_on_target'] = family_data['FamMetrics>>on_target.unique_molecules'] * \
                                                                desired_amplification
    
    family_data['FamMetrics>>molecules_for_optimal_off_target'] = family_data['FamMetrics>>off_target.unique_molecules'] * \
                                                                ( desired_amplification / family_data['FamMetrics>>amplification_bias.on/off'])
    
    family_data['FamMetrics>>total_molecules_for_optimal'] = family_data['FamMetrics>>molecules_for_optimal_on_target'] + family_data['FamMetrics>>molecules_for_optimal_off_target']
    family_data['FamMetrics>>total_GBs_for_optimal'] = family_data['FamMetrics>>total_molecules_for_optimal'] * 300 / 1000000000

    family_data["FamMetrics>>on_target.duplex_conversion_rate"] = family_data['FamMetrics>>molecules_for_optimal_on_target'] / family_data['FamMetrics>>on_target.duplex.dscs']
    family_data["FamMetrics>>off_target.duplex_conversion_rate"] = family_data['FamMetrics>>molecules_for_optimal_off_target'] / family_data['FamMetrics>>off_target.duplex.dscs']
    family_data["FamMetrics>>all.duplex_conversion_rate"] = family_data['FamMetrics>>total_molecules_for_optimal'] / family_data['FamMetrics>>all.duplex.dscs']

    family_data['FamMetrics>>all.GBs'] = family_data['FamMetrics>>all.raw_reads'] * 300 / 1000000000
    family_data['FamMetrics>>missing_GBs_for_optimal'] = family_data['FamMetrics>>total_GBs_for_optimal'] - family_data['FamMetrics>>all.GBs']
    family_data['FamMetrics>>missing_GBs_for_optimal'] = family_data['FamMetrics>>missing_GBs_for_optimal'].apply(lambda x: 0 if x < 0 else x)

    family_data['FamMetrics>>missing_raw_reads_for_optimal'] = family_data['FamMetrics>>total_molecules_for_optimal'] - family_data['FamMetrics>>all.raw_reads']
    family_data['FamMetrics>>missing_raw_reads_for_optimal'] = family_data['FamMetrics>>missing_raw_reads_for_optimal'].apply(lambda x: 0 if x < 0 else x)

    return family_data


# Compute metrics from family metrics

def compute_metrics_from_familymetrics(on_target_file, all_file, on_target_curve_file, all_curve_file):
    data_met = fam_metrics2df(on_target_file)
    data_met_ont = fam_metrics2df(all_file, prefix ='on_target')

    data_curv = fam_curves2df(on_target_curve_file)
    data_curv_ont = fam_curves2df(all_curve_file, prefix = 'on_target')

    data = data_met.merge(data_curv, on='sample')
    data_ont = data_met_ont.merge(data_curv_ont, on='sample')
    
    data_family = data.merge(data_ont, on='sample')
    
    # print(data_family.columns.tolist())

    off_targets_families = compute_off_targets_metrics(data_family)
    data_family_joined = data_family.merge(off_targets_families, on='sample')
    
    return compute_desired_sequencing_metrics(data_family_joined)



@click.command()
@click.option('--on_target_metrics', '-t', required=True, type=click.Path(exists=True), help='Path to the on-target metrics file.')
@click.option('--all_metrics', '-a', required=True, type=click.Path(exists=True), help='Path to the all metrics file.')
@click.option('--on_target_curve', '-c', required=True, type=click.Path(exists=True), help='Path to the on-target curve file.')
@click.option('--all_curve', '-C', required=True, type=click.Path(exists=True), help='Path to the all curve file.')
@click.option('--output-prefix', '-p', required=True, type=click.Path(), help='Prefix for the output files.')
def main(on_target_metrics, all_metrics, on_target_curve, all_curve, output_prefix):
    """
    Main function to process the input file(s) and generate plots.
    Accepts one or more input files for aggregation.
    """

    # Convert tuple to list, or handle single file
    data_family_complete = compute_metrics_from_familymetrics(on_target_metrics, all_metrics, on_target_curve, all_curve)
    
    # store different pieces separately
    # FamMetrics>>off_target
    data_family_complete[["sample"] + [x for x in data_family_complete.columns if x.startswith("FamCurves>>")]
                         ].to_csv(output_prefix + ".FamCurves.tsv", sep='\t', index=False)
    
    data_family_complete[["sample"] + [x for x in data_family_complete.columns if not x.startswith("FamCurves>>")]
                             ].to_csv(output_prefix + ".tsv", sep='\t', index=False)

if __name__ == '__main__':
    main()
