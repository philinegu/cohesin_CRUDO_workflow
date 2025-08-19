#!/usr/bin/env python3

"""
05b_CRUDO_RAD21_analysis.py
Author: Philine Guckelberger
Date: 2025/07/30

Description:
    Reads CRUDO enhancer data and corresponding BED signal files (in this case RAD21 ChIP-seq outputs from _05a_PrepareBeds_CRUDO_RAD21_comp.sh),
    merges them, computes log2 fold-change and % change between treated/untreated,
    and visualizes paired enhancer & CTCF signals with boxplots.

Inputs:
    --CSV file containing CRUDO enhancers (e.g., 04b_CRUDO_FF_enhancers.csv)

Outputs:
    - A CSV with merged signal data
    - A boxplot visualization of paired enhancer & CTCF signals

Usage:
    python scripts/_05b_CRUDO_RAD21_analysis.py \
        --enhancers resources/CRUDO_FF_enhancers.csv \
        --output_directory_csv path/to/output/directory/ \
        --output_directory_plot path/to/output/directory/plots/
"""


# ------------------ Imports ------------------
import argparse
import os
import pandas as pd
import numpy as np
import pybedtools
from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ CONFIGURATION ------------------

# Hardcoded BED files, paths specified in _05a_preparing_beds.sh
element_bed = {
    "path/to/output/directory/CRUDO_hg19_enhancers_noTSS_signal_GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bed": "RAD21.signal.noAux",
    "path/to/output/directory/CRUDO_hg19_enhancers_noTSS_signal_GSM2809610_Rao-2017-CHIP002-RAD21-treated.bed": "RAD21.signal.Aux"
}
ctcf_bed = {
    "path/to/output/directory/CRUDO_CTCF_within_target_gene_loci_hg19_signal_GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bed": "RAD21.signal.noAux",
    "path/to/output/directory/CRUDO_CTCF_within_target_gene_loci_hg19_signal_GSM2809610_Rao-2017-CHIP002-RAD21-treated.bed": "RAD21.signal.Aux"
}

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Merge enhancer/CTCF BED signals and plot paired boxplots.")
    parser.add_argument("--enhancers", required=True, help="CRUDO enhancer csv. Output from scripts/_04b_CRUDO_FF_enhancer_classification.py")
    parser.add_argument("--output_directory_csv", required=True, help="Output diretory for CSVs.")
    parser.add_argument("--output_directory_plot", required=True, help="Output diretory for plots.")
    return parser.parse_args()


def preprocess_signal(file_path, sample_name):
    df = pybedtools.BedTool(file_path).to_dataframe(names=['name','len1','len2','sum','mean','mean2'])
    return df[['name','mean']].rename(columns={'mean': sample_name})


def merge_signals(bed_dict):
    signal_dfs = [preprocess_signal(f, s) for f, s in bed_dict.items()]
    return reduce(lambda l, r: pd.merge(l, r, on='name', how='outer'), signal_dfs)

def paired_boxplot(final_df, output_plot, file_name=None):
    # Prepare for plotting
    data = pd.DataFrame({
        'name': final_df['TargetGene_Category'],
        'noAux': final_df['RAD21.signal.noAux'],
        'plusAux': final_df['RAD21.signal.Aux']
    })
    data_melted = data.melt(id_vars=['name'], var_name='Condition', value_name='Value')
    unique_names = data['name'].unique()
    positions = {name: i for i, name in enumerate(unique_names)}
    # Plot
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(
        x='name', 
        y='Value', 
        hue='Condition', 
        data=data_melted, 
        width=0.6, 
        palette=["#6E7CA0", "#8F3237"], 
        showfliers=False,
        ax=ax
    )
    # Overlay paired points & connecting lines
    for name in unique_names:
        name_data = data[data['name'] == name]
        x_position = positions[name]
        rep1_values = name_data['noAux'].values
        rep2_values = name_data['plusAux'].values
        x_rep1 = [x_position - 0.1] * len(rep1_values)
        x_rep2 = [x_position + 0.1] * len(rep2_values)
        ax.plot(x_rep1, rep1_values, 'o', markersize=3, alpha=0.5, color='#6E7CA0')
        ax.plot(x_rep2, rep2_values, 'o', markersize=3, alpha=0.5, color='#8F3237')
        for x1, x2, y1, y2 in zip(x_rep1, x_rep2, rep1_values, rep2_values):
            ax.plot([x1, x2], [y1, y2], color='gray', linewidth=0.5, alpha=0.6)

    ax.set_xticks([positions[name] for name in unique_names])
    ax.set_xticklabels(unique_names, rotation=30, ha="right")
    ax.set_xlabel('')
    ax.set_ylabel('RAD21 signal')
    plt.axhline(y=0, color='k', linestyle='--', lw=1)
    plt.axhline(y=5, color='k', linestyle='--', lw=1)
    plt.legend(title='Condition')
    plt.tight_layout()
    if file_name:
        path_plot = os.path.join(output_plot, file_name)
        plt.savefig(path_plot, format="pdf", bbox_inches="tight")
        print(f'Saving to: {path_plot}')
    plt.close(fig)

# ------------------ Main ------------------

def main(args):
    # Load enhancer metadata
    enhancers = pd.read_csv(args.enhancers, sep=None, engine="python")

    # Merge enhancer & CTCF signals
    merged_element_df = merge_signals(element_bed)
    merged_ctcf_df = merge_signals(ctcf_bed)

    # Merge with enhancer metadata & compute log2FC
    full_element_df = enhancers.merge(merged_element_df, on='name', how='left')
    full_element_df['signal.Rad21.log2FC'] = np.log2(full_element_df['RAD21.signal.Aux']) - np.log2(full_element_df['RAD21.signal.noAux'])
    full_element_df['%Change.signal.Rad21'] = (full_element_df['RAD21.signal.Aux'] / full_element_df['RAD21.signal.noAux']) - 1

    os.makedirs(args.output_directory_csv, exist_ok=True)
    path_csv = os.path.join(args.output_directory_csv, "CRUDO_FF_enhancers_RAD21.csv")
    full_element_df.to_csv(path_csv, index=False)
    print(f'Saving to: {path_csv}')

    # Annotate categories
    merged_element_df['category'] = 'enhancer'
    merged_ctcf_df['category'] = 'CTCF'
    gene_mapping = {'chr2':'SSFA2','chr11':'CCND1','chr8':'MYC','chr7':'FAM3C','chr12':'KITLG'}
    for df in [merged_element_df, merged_ctcf_df]:
        df['TargetGene'] = df['name'].str.split(':').str[0].map(gene_mapping)

    # Combine enhancer & CTCF data
    final_df = pd.concat([merged_element_df, merged_ctcf_df], ignore_index=True)
    final_df['TargetGene_Category'] = final_df['TargetGene'] + ' - ' + final_df['category']

    # Plot paired box plot (like Extended Data Fig. 11f)
    output_plot = args.output_directory_plot
    os.makedirs(args.output_directory_plot, exist_ok=True)
    paired_boxplot(final_df, output_plot, file_name="Boxplot_RAD21.pdf")

if __name__ == "__main__":
    args = parse_args()
    main(args)