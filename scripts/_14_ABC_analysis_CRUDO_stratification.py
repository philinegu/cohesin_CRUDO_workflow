#!/usr/bin/env python3
"""
14_ABC_analysis_CRUDO_stratification.py
Author: Philine Guckelberger
Date: 2025/08/12

Description:
    Analyze ABC-predicted enhancer-gene paris with and without cohesin 
    (Aux vs. noAux) and compare retention, distance, contact strength, and CTCF 
    binding. Stratify enhancer-gene pairs by CRUDO PRO-seq expression groups and generate 
    summary statistics and plots.


Inputs:
    - A CSV consisting of thresholded ABC predictions in the presence of cohesin (noAux, in this case equivalent to Supplememtary Table 1a)
    - A CSV consisting of thresholded ABC predictions in the absence of cohesin (Aux, in this case equivalent to Supplememtary Table 1b)
    - A PRO-seq TPM txt file (in this case, output from _13_CRUDO_PROseq_GeneExpression.py)
    - A PRO-seq DeSeq2 txt file (in this case, Auxin_vs_Control.RAD21.Genes.DESeq2.txt)
    - A list of housekeeping genes (defined in Fulco et al. 2019 Supplementary Table 2d)
    - A CSV file containing CRUDO enhancers and extracted feature data (output from _05b_CRUDO_RAD21_analysis.py)

Outputs:
    - Rolling average plots of Hi-C contact, ABC score, and enhancer retention (and scatter plots of raw data)
    - Boxplots and barplots comparing enhancer categories across expression groups
    - Per-gene ABC enhancer summary tables
    - Statistics and summary plots on retained vs. lost enhancers, odds ratios, and expression-group distributions


Usage:
    python scripts/_14_ABC_analysis_CRUDO_stratification.py \
        --predictions_noAux  resources/1a_ABC_thresholded_noAux.csv\
        --predictions_Aux  resources/1b_ABC_thresholded_Aux.csv\
        --PRO_TPM resources/CRUDO_Genes_Pro.tpm.txt\
	    --PRO_DeSeq2 resources/Auxin_vs_Control.RAD21.Genes.DESeq2.txt\
        --housekeeping_genes  resources/HousekeepingGenes_Fulco2019_S5d.txt\
	    --CRUDO_enhancers  resources/CRUDO_FF_enhancers_RAD21.csv\
        --output_directory_plot path/to/output/directory/plots/

"""


# ------------------ Imports ------------------
import argparse
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from matplotlib.ticker import FormatStrFormatter
from scipy import stats
from scipy.stats import gaussian_kde
import math
from scipy.stats.contingency import odds_ratio
from scipy.stats import fisher_exact
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Analyze ABC enhancer-gene predictions with/without cohesin and stratify by CRUDO expression")
    parser.add_argument("--predictions_noAux", help="Thresholded ABC predictions in the presence of cohesin (ie SupplementaryTable 1a)")
    parser.add_argument("--predictions_Aux", help="Thresholded ABC predictions in the absence of cohesin (ie SupplementaryTable 1b)")
    parser.add_argument("--PRO_TPM", help="Path to input txt file with PRO-seq TPM data (ie output _13_CRUDO_PROseq_GeneExpression.py)")
    parser.add_argument("--PRO_DeSeq2", help="Path to input txt file with PRO-seq DeSeq2 values")
    parser.add_argument("--housekeeping_genes", help="List of housekeeping genes for filtering (ie Fulco et al. 2019 SupplementaryTable 5d)")
    parser.add_argument("--CRUDO_enhancers", help="CRUDO enhancer csv (output from _05b_CRUDO_RAD21_analyzsis.py)")
    parser.add_argument("--output_directory_plot", help="Directory to save outputs")
    return parser.parse_args()


def add_relevant_columns(df):
    # Scale Hi-C AUX counts based on overall contacts in each conditon
    total_noAux = 1093601688
    total_Aux = 1527291641
    scaling_factor = total_noAux / total_Aux
    df['SCALE.normalized.observed.5Kb.Aux.scaled'] = df['SCALE.normalized.observed.5Kb.Aux'] * scaling_factor
    #Add columns for HiC, ABC, and H3K27ac changes
    # Hi-C changes
    df['normalized.HiC.5Kb.noAux.pseudo'] = df['SCALE.normalized.observed.5Kb.noAux'] + 1
    df['normalized.HiC.5Kb.Aux.pseudo'] = df['SCALE.normalized.observed.5Kb.Aux.scaled'] + 1
    df['normalized.HiC.5Kb.log2FC.(Aux/noAux)'] = np.log2(df['normalized.HiC.5Kb.Aux.pseudo']) - np.log2(df['normalized.HiC.5Kb.noAux.pseudo'])
    hic_FC = df['normalized.HiC.5Kb.Aux.pseudo'] / df['normalized.HiC.5Kb.noAux.pseudo']
    df['%Change.normalized.HiC.5Kb'] = (hic_FC - 1)
    # ABC changes
    ABC_FC = df['ABC.Score.Aux'] / df['ABC.Score.noAux']
    df['ABC.Score.log2FC.(Aux/noAux)'] = np.log2(df['ABC.Score.Aux']) - np.log2(df['ABC.Score.noAux'])
    df['%Change.ABC.Score'] = (ABC_FC - 1)
    # H3K27ac changes
    H3K27ac_FC = df['H3K27ac.RPM.values.Aux'] / df['H3K27ac.RPM.values.noAux']
    df['H3K27ac.log2FC.(Aux/noAux)'] = np.log2(df['H3K27ac.RPM.values.Aux']) - np.log2(df['H3K27ac.RPM.values.noAux'])
    df['%Change.H3K27ac.RPM.values'] = (H3K27ac_FC - 1)
    return df


def filter_predictions(df, df_tpm, df_hk, tpm_threshold=1, remove_promoters=True, remove_hk=True):
    df_tpm_filtered = df_tpm.loc[df_tpm['RAD21_NT_avg'] > tpm_threshold]
    filtered = df[df['TargetGene'].isin(df_tpm_filtered['Gene Symbol'])]
    #remove TSS elements
    if remove_promoters:
        filtered = filtered[filtered['class'] != 'promoter']
    #remove HK genes
    if remove_hk:
        filtered = filtered[~filtered['TargetGene'].isin(df_hk['GeneSymbol'])]
    
    expressed_nonHK= df_tpm_filtered [~df_tpm_filtered ['Gene Symbol'].isin(df_hk['GeneSymbol'])]
    
    return filtered, df_tpm_filtered, expressed_nonHK


def compute_percent_retained(df_present, df_absent):
    enhancers_noAux = df_present['name'].unique()
    enhancers_Aux = df_absent['name'].unique()
    df_present['PercentRetained'] = df_present['name'].apply(
        lambda e: 100 if e in enhancers_Aux else 0
         )
    return df_present


def enhancer_stats(df):
    genes = df['TargetGene'].unique()
    n_enhancers = []
    dist_mean = []
    dist_median = []
    for gene in genes:
        subset = df.loc[df['TargetGene'] == gene, 'DistanceToTSS.Kb']
        n_enhancers.append(len(subset))
        dist_mean.append(np.mean(subset))
        dist_median.append(np.median(subset))
    stats_summary = {
        'num_genes': len(genes),
        'num_pairs': len(df),
        'avg_enhancers_per_gene': np.mean(n_enhancers),
        'avg_mean_distance_per_gene': np.mean(dist_mean),
        'avg_median_distance_per_gene': np.mean(dist_median),
        'total_mean_distance': np.mean(df['DistanceToTSS.Kb']),
        'total_median_distance': np.median(df['DistanceToTSS.Kb']),
        'p90_distance': np.percentile(df['DistanceToTSS.Kb'], 90)
    }
    return stats_summary

def summarize_columns(df, condition, column, label=""):
    subset = df.loc[condition]
    print(f"{label} count:", len(subset))
    print(f"{label} mean {column} ", subset[column].mean())
    print("")


def compute_rolling_metrics(df, sort_by='DistanceToTSS.Kb'):
    # Sort
    df_sorted = df.sort_values(by=[sort_by], ascending=True).copy()
    # Define column groups: col -> (window, min_periods)
    col_settings = {
        'SCALE.normalized.observed.5Kb.noAux': (1000, 100),
        'SCALE.normalized.observed.5Kb.Aux.scaled': (1000, 100),
        'normalized.HiC.5Kb.log2FC.(Aux/noAux)': (1000, 100),
        'ABC.Score.noAux': (1000, 100),
        'ABC.Score.Aux': (1000, 100),
        'ABC.Score.log2FC.(Aux/noAux)': (1000, 100),
        'PercentRetained': (1000, 100),
    }
    # Compute rolling stats for each
    for col, (window, min_p) in col_settings.items():
        roll_mean = df_sorted[col].rolling(window, min_periods=min_p).mean()
        roll_std = df_sorted[col].rolling(window, min_periods=min_p).std()

        df_sorted[f'Rolling_{col}'] = roll_mean
        df_sorted[f'Rolling_{col}.stdev'] = roll_std
        df_sorted[f'Rolling_{col}.ci95'] = 1.96 * (roll_std / np.sqrt(window))
     #Filter to pnlu plot enhancers between 1-1000
    df_sorted = df_sorted.loc[(df_sorted['DistanceToTSS.Kb']>=1) & (df_sorted['DistanceToTSS.Kb']<=1000)]
    
    return df_sorted

def plot_rolling_with_ci(
    data,
    x_col,
    y_noAux_col,
    y_Aux_col,
    ci_noAux_col,
    ci_aux_col,
    color_noAux='#6E7CA0',
    color_aux='#8F3237',
    rasterize=False,
    save_path=None,
    ylim=None
    ):
    fig, ax = plt.subplots(figsize=(4.25,1.25))
    # Plot main lines
    l1 = sns.lineplot(data=data, x=x_col, y=y_noAux_col, color=color_noAux, errorbar=None, ax=ax)
    l2 = sns.lineplot(data=data, x=x_col, y=y_Aux_col, color=color_aux, errorbar=None, ax=ax)
    # Extract values
    x = data[x_col]
    y_noAux = data[y_noAux_col]
    y_Aux = data[y_Aux_col]
    # Plot error bands
    fill1 = ax.fill_between(x, y_noAux - data[ci_noAux_col], y_noAux + data[ci_noAux_col], color=color_noAux, alpha=0.2)
    fill2 = ax.fill_between(x, y_Aux - data[ci_aux_col], y_Aux + data[ci_aux_col], color=color_aux, alpha=0.2)
    # Optional rasterization
    if rasterize:
        for line in ax.lines:
            line.set_rasterized(True)
        fill1.set_rasterized(True)
        fill2.set_rasterized(True)
    # Log-log scaling and formatting
    ax.set_xlim(1, 1000)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Save figure if requested
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    return fig, ax

def plot_single_rolling_with_ci(
    data,
    x_col,
    y_col,
    ci_col,
    color='k',
    rasterize=False,
    save_path=None,
    ylim=None
    ):
    # Define plot
    fig, ax = plt.subplots(figsize=(4.25,1.25))
    # Plot main line
    line = sns.lineplot(data=data, x=x_col, y=y_col, color=color, errorbar=None, ax=ax)
    # CI fill
    y = data[y_col]
    yerr0 = y - data[ci_col]
    yerr1 = y + data[ci_col]
    fill = ax.fill_between(data[x_col], yerr0, yerr1, color=color, alpha=0.2)
    # Rasterization option
    if rasterize:
        for l in ax.lines:
            l.set_rasterized(True)
        fill.set_rasterized(True)
    # Axes formatting
    ax.set_xlim(1, 1000)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if ylim:
        ax.set_ylim(*ylim)
    # Add dashed line at y 0 (highlights a fold change of 1)
    ax.axhline(y=0, xmin=-1, xmax=1, color='black', lw=1, linestyle='dashed')
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Save
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    return fig, ax

def plot_rasterized_scatter_dual(
    data,
    x_col,
    y_noAux_col,
    y_Aux_col,
    color_left='#6E7CA0',
    color_right='#8F3237',
    rasterize=True,
    save_path=None,
    ylim=None,
    ):
    # Define plots
    fig, axes = plt.subplots(1, 2, figsize=(8.5, 1.75), sharey=True)
    # Left panel scatter
    axes[0].scatter(data[x_col], data[y_noAux_col], color=color_left, alpha=0.2, s=8)
    # Right panel scatter
    axes[1].scatter(data[x_col], data[y_Aux_col], color=color_right, alpha=0.2, s=8)
    # Rasterize scatter points only
    if rasterize:
        for ax in axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
    # Shared formatting
    for ax in axes:
        ax.set_xlim(1, 1000)
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        if ylim:
            ax.set_ylim(*ylim)
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        sns.despine(ax=ax)
    #Save plot
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format='pdf', bbox_inches='tight')
    plt.close(fig)

    return fig, axes


def plot_rasterized_scatter(
    data,
    x_col,
    y_col,
    color='k',
    rasterize=True,
    save_path=None,
    ylim=None
    ):
    # Define plot
    fig, ax = plt.subplots(figsize=(4.25,1.25))
    # Plot scatter
    ax = sns.scatterplot(data=data, x=x_col, y=y_col, color=color, s=8, ax=ax)
    # Rasterize scatter points only
    if rasterize: 
        for coll in ax.collections:
            coll.set_rasterized(True)
    # Axes formatting
    ax.set_xlim(1, 1000)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if ylim:
        ax.set_ylim(*ylim)
    # Add dashed line at y 0 (highlights a fold change of 1)
    ax.axhline(y=0, xmin=-1, xmax=1, color='black', lw=1, linestyle='dashed')
    # Hide right and top spines
    ax.spines[['right', 'top']].set_visible(False)
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
        
    plt.close(fig)
    return fig, ax


def plot_histogram_rasterized(
    data,
    col,
    color='k',
    rasterize=True,
    save_path=None,
    ylim=None
    ):
   # Define plot
    fig, ax = plt.subplots(figsize=(4.25,1.25))
    # Plot histogramm
    ax = sns.histplot(data[col], color=color, log_scale=True, element="step", fill=False)
    # Rasterize histogram patches (outlines)
    if rasterize:
        for patch in ax.patches:
            patch.set_rasterized(True)
    # Axes formatting
    ax.set_xlim(1, 1000)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    #Save plot
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)

    return fig, ax


def plot_hic_contact_rolling(df, output_dir, filename=None):
    # plots Hi-C contact (cohesin present and absent) as a rolling circle average (like Fig. 1a)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None

    fig, ax = plot_rolling_with_ci(
        df,
        x_col="DistanceToTSS.Kb",
        y_noAux_col="Rolling_SCALE.normalized.observed.5Kb.noAux",
        y_Aux_col="Rolling_SCALE.normalized.observed.5Kb.Aux.scaled",
        ci_noAux_col="Rolling_SCALE.normalized.observed.5Kb.noAux.ci95",
        ci_aux_col="Rolling_SCALE.normalized.observed.5Kb.Aux.scaled.ci95",
        rasterize=True,
        ylim=(1,1000),
        save_path=save_path
        )
    # Calculate power law values for 'DistanceToTSS.Kb'
    df['powerlaw_distance'] = 500 * df['DistanceToTSS.Kb'] ** -1
    # Plot power law values
    sns.lineplot(data=df, x="DistanceToTSS.Kb", y="powerlaw_distance", color='black', errorbar=None, ax=ax)
    plt.close(fig)
    return fig, ax


def plot_abc_score_rolling(df, output_dir, filename=None):
# plots ABC score (cohesin present and absent) as a rolling circle average (like Fig. 1c)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:   
        save_path = None
    fig, ax = plot_rolling_with_ci(
        df,
        x_col="DistanceToTSS.Kb",
        y_noAux_col="Rolling_ABC.Score.noAux",
        y_Aux_col="Rolling_ABC.Score.Aux",
        ci_noAux_col="Rolling_ABC.Score.noAux.ci95",
        ci_aux_col="Rolling_ABC.Score.Aux.ci95",
        rasterize=True,
        save_path=save_path
         )
    ax.axhline(y=0.025693, xmin=-1, xmax=1, color='black', lw=1, label=False, linestyle='dashed')
    ax.yaxis.set_minor_formatter(FormatStrFormatter("%.3f"))
    plt.close(fig)
    return fig, ax


def plot_log2fc_hic_contact_rolling(df, output_dir, filename=None):
    # plots change in Hi-C contact as a rolling average (like Fig. 1b)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax= plot_single_rolling_with_ci(
        data=df,
        x_col="DistanceToTSS.Kb",
        y_col="Rolling_normalized.HiC.5Kb.log2FC.(Aux/noAux)",
        ci_col="Rolling_normalized.HiC.5Kb.log2FC.(Aux/noAux).ci95",
        color='k',
        rasterize=True,
        ylim=(-2.5,0),
        save_path=save_path
         )
    plt.close(fig)
    return fig, ax


def plot_log2fc_abc_score_rolling(df, output_dir, filename=None):
    # plots change in Hi-C contact as a rolling average (like Fig. 1d)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax= plot_single_rolling_with_ci(
        data=df,
        x_col="DistanceToTSS.Kb",
        y_col="Rolling_ABC.Score.log2FC.(Aux/noAux)",
        ci_col="Rolling_ABC.Score.log2FC.(Aux/noAux).ci95",
        color='k',
        rasterize=True,
        ylim=(-2.5,0),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax


def plot_log2fc_percent_retained_rolling(df, output_dir, filename=None):
    # plots percent enhancers retained in the absence of cohesin as a rolling average (like Fig. 1g)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax= plot_single_rolling_with_ci(
        data=df,
        x_col="DistanceToTSS.Kb",
        y_col="Rolling_PercentRetained",
        ci_col="Rolling_PercentRetained.ci95",
        color='k',
        rasterize=True,
        ylim=(0,100),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax


def plot_hic_contact_scatter(df, output_dir, filename=None):
    # plots Hi-C contact (cohesin present and absent) scatter dual plot (like Extended Data Fig. 1a)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, axes = plot_rasterized_scatter_dual(
        df,
        x_col='DistanceToTSS.Kb',
        y_noAux_col="SCALE.normalized.observed.5Kb.noAux",
        y_Aux_col="SCALE.normalized.observed.5Kb.Aux.scaled",
        rasterize=True,
        ylim=(0.1,1000),
        save_path=save_path
        )
    plt.close(fig)
    return fig, axes

def plot_abc_score_scatter(df, output_dir, filename=None):
    # plots ABC score scatter (cohesin present and absent) dual plot (like Extended Data Fig. 1i)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, axes = plot_rasterized_scatter_dual(
        df,
        x_col="DistanceToTSS.Kb",
        y_noAux_col="ABC.Score.noAux",
        y_Aux_col="ABC.Score.Aux",
        rasterize=True,
        save_path=save_path
        )
    for ax in axes:
        ax.axhline(y=0.025693, xmin=0, xmax=1, color='black', lw=1, label=False, linestyle='dashed')
    plt.close(fig)
    return fig, axes


def plot_log2fc_hic_contact_scatter(df, output_dir, filename=None):
    # plots change in Hi-C contact scatter (like Extended Data Fig. 1b)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax = plot_rasterized_scatter(
        df,
        x_col='DistanceToTSS.Kb',
        y_col='normalized.HiC.5Kb.log2FC.(Aux/noAux)',
        rasterize=True,
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax


def plot_log2fc_abc_score_scatter(df, output_dir, filename=None):
    # plots change in ABC score scatter (like Extended Data Fig. 1j)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax = plot_rasterized_scatter(
        df,
        x_col='DistanceToTSS.Kb',
        y_col='ABC.Score.log2FC.(Aux/noAux)',
        rasterize=True,
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax

def plot_distance_histogram(df, output_dir, filename=None):
    # plots histogram of DistanceToTSS (like Extended Data Fig. 1f)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    fig, ax = plot_histogram_rasterized(
        df,
        "DistanceToTSS.Kb",
        color='k',
        rasterize=True,
        save_path = save_path,
        ylim=(0,1000)
        )
    yticks = np.arange(0, 1001, 200)
    ax.set_yticks(yticks)
    plt.close(fig)
    return fig, ax

def plot_abc_enhancers_cdf(
    df_filtered, 
    df_aux_filtered, 
    gene_col='TargetGene', 
    dist_col='DistanceToTSS.Kb', 
    label_present='Cohesin present', 
    label_absent='Cohesin absent',
    color_present='#6E7CA0',
    color_absent='#8F3237',
    save_path_enhancers=None,
    save_path_distance=None
    ):
    # Compute number of enhancers per gene for both conditions
    n_enhancers = df_filtered.groupby(gene_col).size().values
    n_enhancers_aux = df_aux_filtered.groupby(gene_col).size().values
    # Plot CDF of number of enhancers per gene
    fig, ax = plt.subplots(figsize=(5,5))
    sns.ecdfplot(n_enhancers, label=label_present, color=color_present, lw=4)
    sns.ecdfplot(n_enhancers_aux, label=label_absent, color=color_absent, lw=4)
    # Set axes
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel("n_enhancers")
    ax.set_ylabel("CDF")
    ax.legend()
    print("Mann-Whitney U test for number of enhancers per gene:", 
          stats.mannwhitneyu(n_enhancers, n_enhancers_aux))
    print(f"Number of genes (cohesin present): {len(n_enhancers)}")
    print(f"Number of genes (cohesin absent): {len(n_enhancers_aux)}")
    if save_path_enhancers:
        plt.savefig(save_path_enhancers, format='pdf', bbox_inches='tight')
    plt.close(fig)
    # Plot CDF of enhancer distance to TSS
    fig, ax = plt.subplots(figsize=(5,5))
    sns.ecdfplot(df_filtered[dist_col], label=label_present, color=color_present, lw=4)
    sns.ecdfplot(df_aux_filtered[dist_col], label=label_absent, color=color_absent, lw=4)
    # Set aces
    ax.set_xlim(left=1)
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_xlabel(dist_col)
    ax.set_ylabel("CDF")
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()
    print("Mann-Whitney U test for enhancer distances:", 
          stats.mannwhitneyu(df_filtered[dist_col], df_aux_filtered[dist_col]))
    if save_path_distance:
        plt.savefig(save_path_distance, format='pdf', bbox_inches='tight')
    plt.close(fig)



def plot_boxplot_from_dict(data_dict, 
                          ylabel="Value", 
                          palette="rocket", 
                          yscale=None,
                          ylim=None, 
                          hline=0,  
                          figsize=(5,5), 
                          save_path=None):
    # Combine data into a DataFrame
    combined_df = []
    for label, data in data_dict.items():
        temp_df = pd.DataFrame({'value': data})
        temp_df['data'] = label
        combined_df.append(temp_df)
    combined_df = pd.concat(combined_df).reset_index(drop=True)
    #Define plot
    fig, ax = plt.subplots(figsize=(5,5))
    # Make boxplot
    ax = sns.boxplot(x="data", y="value", data=combined_df, palette=palette, fliersize=0.01, hue = "data", legend = False, 
                     medianprops={"color": "black"})
    #Plot dashed line
    if hline is not None:
        ax.axhline(y=hline, xmin=0, xmax=1, color='black', lw=1, linestyle="dashed")
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Set axes
    if yscale:
        ax.set_yscale(yscale)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if ylim:
        ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    # Save plot
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches='tight')

    plt.close(fig)
    return fig, ax


def plot_basic_boxplot(df, group_col, value_col, palette="rocket", 
                       hline=0, yscale=None,  ylim=None,  order=None, save_path=None):
    #Define plot
    fig, ax = plt.subplots(figsize=(5,5))
    # Make boxplot
    ax = sns.boxplot(
        x=group_col,
        y=value_col,
        data=df,
        palette=palette,
         hue = group_col, legend = False, 
        showfliers=False,
        medianprops={"color": "black"},
        order=order
        )
    #Plot dashed line
    if hline is not None:
        ax.axhline(y=hline, xmin=0, xmax=1, color='black', lw=1, linestyle="dashed")
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Set axes
    if yscale:
        ax.set_yscale(yscale)
        ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if ylim:
        ax.set_ylim(*ylim)
    # Save plot
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches='tight')
    plt.close(fig)
    return fig, ax


def plot_h3k27_hic_boxplot(df, output_dir, filename=None):
    # Plots boxplot of log2FC in H3K27ac and hi-c contact across ABC enhancers (like Extended Data Figure 1d)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    # Create data dict for plotting
    dict_h3k27_hic = {
        "h3k27ac": df['H3K27ac.log2FC.(Aux/noAux)'].dropna(),
        "hic": df['normalized.HiC.5Kb.log2FC.(Aux/noAux)'].dropna()
        }   
    # Generate plot
    fig, ax = plot_boxplot_from_dict(
        dict_h3k27_hic,
        ylabel="log2 Fold Change",
        ylim=(-5, 5),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax

def plot_retained_hic_boxplot(df, output_dir, filename=None):
    # Plot boxplot of log2FC in hi-c contact for retained vs lost ABC enhancers (like Extended Data Figure 1e)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    #Create data dict for plotting
    dict_retained_hic = {
        "retained": df.loc[df['PercentRetained']==100]['normalized.HiC.5Kb.log2FC.(Aux/noAux)'],
        "lost": df.loc[df['PercentRetained']==0]['normalized.HiC.5Kb.log2FC.(Aux/noAux)']
        }   
    # Generate plot
    fig, ax =plot_boxplot_from_dict(
        dict_retained_hic,
        ylabel="normalized.HiC.5Kb.log2FC.(Aux/noAux))",
        palette=['#6D5D6E','#F4EEE0'],
        ylim=(-5, 2.5),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax

def plot_retained_dist_boxplot(df, output_dir, filename=None):
    # Plot boxplot of distance to TSS for retained vs lost ABC enhancers (like Extended Data Figure 1f)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    #Create data dict for plotting
    dict_retained_distance = {
        "retained": df.loc[df['PercentRetained']==100]['DistanceToTSS.Kb'],
        "lost": df.loc[df['PercentRetained']==0]['DistanceToTSS.Kb']
            }   
    # Generate plot
    fig, ax =plot_boxplot_from_dict(
        dict_retained_distance,
        ylabel="DistanceToTSS.Kb",
        palette=['#6D5D6E','#F4EEE0'],
        yscale=('log'),
        ylim=(0.5, 1000),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax


def plot_cohesin_distance_boxplot(df_noAux, df_Aux, output_dir, filename=None):
    # Plot boxplot of distance to TSS for ABC predictions in cohesin present and absent conditions (like Extended Data Figure 1h)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    #Create data dict for plotting
    dict_distance = {
        "cohesin present": df_noAux['DistanceToTSS.Kb'],
        "cohesin absent": df_Aux['DistanceToTSS.Kb']
        }   
    # Generate plot 
    fig, ax =plot_boxplot_from_dict(
        dict_distance,
        ylabel="DistanceToTSS.Kb",
        palette=['#6E7CA0','#8F3237'],
        yscale=('log'),
        ylim=(0.5, 1000),
        save_path=save_path
        )
    plt.close(fig)
    return fig, ax

def plot_CTCF_ABC_boxplot(df, output_dir, filename=None):
    # Plot boxplot of change in ABC score for ABC predictions close to a CTCF binding site o rnot (like Extended Data Figure 1i)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    # Generate plot
    fig, ax = plot_basic_boxplot(df, 'CTCFwithin5Kb', "ABC.Score.log2FC.(Aux/noAux)",
                   palette=['pink','#84A7A1'],  ylim=(-2.5, 2.5),  
                   order=[True, False], hline=0, save_path=save_path)
    plt.close(fig)
    return fig, ax


def plot_CTCF_hic_boxplot(df, output_dir, filename=None):
    # Plot boxplot of change in hi-c contact for ABC predictions close to a CTCF binding site o rnot (like Extended Data Figure 1k)
    if filename:
        save_path = os.path.join(output_dir, filename + ".pdf")
    else:
        save_path = None
    # Generate plot
    fig, ax = plot_basic_boxplot(df, 'CTCFwithin5Kb', "normalized.HiC.5Kb.log2FC.(Aux/noAux)",
                   palette=['pink','#84A7A1'],  ylim=(-3, 3),  
                   order=[True, False], hline=0, save_path=save_path)
    plt.close(fig)
    return fig, ax



def plot_density_correlation(df, x_col, y_col, xlim=None, ylim=None, save_path=None):
    # Drop NA values from both columns
    data = df[[x_col, y_col]].dropna()
    # Calculate density
    xy = np.vstack([data[x_col], data[y_col]])
    z = gaussian_kde(xy)(xy)
    # Plot
    fig, ax = plt.subplots(figsize=(5, 5))
    sc = ax.scatter(data[x_col], data[y_col], c=z, s=0.5, cmap=sns.color_palette("rocket", as_cmap=True))
    # Set axes
    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)
    # Draw horizontal dashed line at y=0 if 0 in ylim range
    if ylim != "" and ylim[0] <= 0 <= ylim[1]:
        ax.axhline(0, linestyle='--', color='gray')
    # Draw vertical dashed line at x=0 if 0 in xlim range
    if xlim != "" and xlim[0] <= 0 <= xlim[1]:
        ax.axvline(0, linestyle='--', color='gray')
    # Draw diagonal only if both limits are provided and they are exactly the same
    if xlim != "" and ylim != "" and xlim == ylim:
        plt.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], linestyle='--', color='gray')
    # Clean plot
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.spines[['right', 'top']].set_visible(False)
    # Pearson correlation
    r, pval = stats.pearsonr(data[x_col], data[y_col])
    print(f"log2FC HiC vs log2FC ABC Pearson R = {r:.2f}, p = {pval:.2e}")
    # Save plot
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
    return fig, ax



def plot_ctcf_distal_proximal_retention(df, distance_col='DistanceToTSS.Kb', ctcf_col='CTCFwithin5Kb',
                                           distal_threshold=50, save_path=None):
    # Split distal and proximal
    df_distal = df.loc[df[distance_col] >= distal_threshold]
    df_prox = df.loc[df[distance_col] < distal_threshold]
    # Split by CTCF presence
    df_plus_distal = df_distal.loc[df_distal[ctcf_col] == True]
    df_no_distal = df_distal.loc[df_distal[ctcf_col] == False]
    df_plus_prox = df_prox.loc[df_prox[ctcf_col] == True]
    df_no_prox = df_prox.loc[df_prox[ctcf_col] == False]
    # define the function to compute fractions
    def calc_fractions(sub_df):
        total = len(sub_df)
        retained = len(sub_df.loc[sub_df['PercentRetained'] == 100])
        lost = len(sub_df.loc[sub_df['PercentRetained'] == 0])
        retained_fraction = retained / total if total > 0 else 0
        lost_fraction = lost / total if total > 0 else 0
        return total, retained, lost, retained_fraction, lost_fraction
    # call function to compute fractions
    total_plus_distal, retained_plus_distal, lost_plus_distal, retained_fraction_plus_distal, lost_fraction_plus_distal = calc_fractions(df_plus_distal)
    total_no_distal, retained_no_distal, lost_no_distal, retained_fraction_no_distal, lost_fraction_no_distal = calc_fractions(df_no_distal)
    total_plus_prox, retained_plus_prox, lost_plus_prox, retained_fraction_plus_prox, lost_fraction_plus_prox = calc_fractions(df_plus_prox)
    total_no_prox, retained_no_prox, lost_no_prox, retained_fraction_no_prox, lost_fraction_no_prox = calc_fractions(df_no_prox)
    # Define labels and fractions to plot:
    labels = ['Dist +CTCF', 'Dist -CTCF', 'Prox +CTCF', 'Prox -CTCF']
    retained_fractions = [retained_fraction_plus_distal, retained_fraction_no_distal,
                             retained_fraction_plus_prox, retained_fraction_no_prox]
    lost_fractions = [lost_fraction_plus_distal, lost_fraction_no_distal,
                         lost_fraction_plus_prox, lost_fraction_no_prox]
    # Define plot
    fig, ax = plt.subplots(figsize=(5,5))
    # Geneerate plot
    ax.bar(labels, retained_fractions, label='Retained', color='#6D5D6E')
    ax.bar(labels, lost_fractions, bottom=retained_fractions, label='Lost', color='#F4EEE0')
    # Clean plot
    ax.set_ylabel('Fraction')
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylim(0, 1)
    ax.legend()
    # Save plot
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
        plt.close(fig)
    return fig, ax


def filter_proseq_and_hic(df_PROseq, df_TPM, df_filtered, TPM_threshold=5, log2fc_dr_threshold=-(np.log2(1.8)),
                           log2fc_non_threshold=-(np.log2(1.4)), padj_threshold=0.05, hic_threshold=5):
    # re-set index
    df_PROseq = df_PROseq.reset_index()
    # Filter genes with RAD21 expression above threshold
    genes_of_interest = df_TPM.loc[df_TPM['RAD21_NT_avg'] > TPM_threshold, 'Gene Symbol']
    df_PROseq_filtered = df_PROseq[df_PROseq['index'].isin(genes_of_interest)]
    # Strongly downregulated
    df_PROseq_filtered_dr = df_PROseq_filtered[
        (df_PROseq_filtered['log2FoldChange'] <= log2fc_dr_threshold) & (df_PROseq_filtered['padj'] <= padj_threshold)
    ]
    # Moderately downregulated
    df_PROseq_filtered_mid = df_PROseq_filtered[
        (df_PROseq_filtered['log2FoldChange'] > log2fc_dr_threshold) &
        (df_PROseq_filtered['log2FoldChange'] <= log2fc_non_threshold) &
        (df_PROseq_filtered['padj'] <= padj_threshold)
    ]
    list_downregulated = df_PROseq_filtered_dr['index'].unique()
    # Filter df_filtered for high HiC signal
    df_filtered_hic = df_filtered[df_filtered['SCALE.normalized.observed.5Kb.noAux'] >= hic_threshold]
    # Define filtered subsets in df_filtered_hic based on PRO-seq downregulation categories
    df_filtered_non = df_filtered_hic[
        (~df_filtered_hic['TargetGene'].isin(df_PROseq_filtered_dr['index'])) &
        (~df_filtered_hic['TargetGene'].isin(df_PROseq_filtered_mid['index'])) &
        (df_filtered_hic['TargetGene'].isin(df_PROseq_filtered['index']))
    ]
    df_filtered_mid = df_filtered_hic[df_filtered_hic['TargetGene'].isin(df_PROseq_filtered_mid['index'])]
    df_filtered_DR = df_filtered_hic[df_filtered_hic['TargetGene'].isin(df_PROseq_filtered_dr['index'])]
    # Print number of genes in each category
    non_fc_value = 2 ** abs(log2fc_non_threshold)
    dr_fc_value = 2 ** abs(log2fc_dr_threshold)
    print(f"Non-Downregulated genes: {len(df_filtered_non['TargetGene'].unique())} (threshold: FC < {non_fc_value:.1f})")
    print(f"Moderately downregulated genes: {len(df_filtered_mid['TargetGene'].unique())} (threshold: FC {non_fc_value:.1f} - {dr_fc_value:.1f})")
    print(f"Strongly downregulated genes: {len(df_filtered_DR['TargetGene'].unique())} (threshold: FC > {dr_fc_value:.1f})")
    print("")
    # Print number of enhancers  with >50% decreace in hic
    count_non = len(df_filtered_non.loc[df_filtered_non['%Change.normalized.HiC.5Kb'] <= -0.5])
    count_mid = len(df_filtered_mid.loc[df_filtered_mid['%Change.normalized.HiC.5Kb'] <= -0.5])
    count_DR = len(df_filtered_DR.loc[df_filtered_DR['%Change.normalized.HiC.5Kb'] <= -0.5])
    print(f"# of enhancers for non-downregulated genes with >50% decrease in hi-c contacts: {count_non} (/ {len(df_filtered_non)}  total)")
    print(f"# of enhancers for moderately downregulated genes with >50% decrease in hi-c contacts: {count_mid} (/ {len(df_filtered_mid)}  total)")
    print(f"# of enhancers for strongly downregulated genes with >50% decrease in hi-c contacts: {count_DR} (/ {len(df_filtered_DR)}  total)")
    print(f"Sum of enhancers with >50% decrease in hi-c contacts: {count_non + count_mid + count_DR}")
    print("")
    return df_filtered_hic, df_filtered_non, df_filtered_mid, df_filtered_DR



def plot_cdf_comparison(df_filtered_non, df_filtered_mid, df_filtered_dr, column, 
                        xlim=None, log_scale=False, threshold=None, vline=None,
                     colors=None, labels=None, save_path=None):
    # Set colors and labels
    colors = colors or ['#646363', 'pink', '#F8B500']
    labels = labels or ["Cohesin-independent genes", "mid genes", "Cohesin-dependent genes"]
    # Plot
    fig, ax = plt.subplots(figsize=(5,5))
    sns.ecdfplot(df_filtered_non[column], color=colors[0], label=labels[0], lw=4, ax=ax)
    sns.ecdfplot(df_filtered_mid[column], color=colors[1], label=labels[1], lw=4, ax=ax)
    sns.ecdfplot(df_filtered_dr[column], color=colors[2], label=labels[2], lw=4, ax=ax)
    # Set axes
    if xlim:
        ax.set_xlim(*xlim)
    if log_scale:
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    if vline is not None:
        ax.axvline(x=vline, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1)
    ax.set_xlabel(column)
    ax.set_ylabel("CDF")
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()
    # Save ploy
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    # Summaries
    n_non = len(df_filtered_non)
    n_dr = len(df_filtered_dr)
    print(f"Total E-P pairs non: {n_non}")
    print(f"Total E-P pairs dr: {n_dr}")

    n_non_thresh = len(df_filtered_non.loc[df_filtered_non[column] < threshold])
    n_dr_thresh = len(df_filtered_dr.loc[df_filtered_dr[column] < threshold])
    print(f"Fraction < {threshold} {column} non: {n_non_thresh / n_non:.3f}")
    print(f"Fraction < {threshold} {column}  kb dr: {n_dr_thresh / n_dr:.3f}")
    print(f"75th percentile non: {np.percentile(df_filtered_non[column], 75):.2f}")
    print(f"75th percentile dr: {np.percentile(df_filtered_dr[column], 75):.2f}")
    # Statistical tests between non and dr groups
    ks_res = stats.ks_2samp(df_filtered_non[column], df_filtered_dr[column], alternative='two-sided', mode='exact')
    ks_p = ks_res.pvalue
    ks_base, ks_exp = f"{ks_p:.2e}".split("e")
    mw_res = stats.mannwhitneyu(df_filtered_non[column], df_filtered_dr[column])
    mw_p = mw_res.pvalue
    mw_base, mw_exp = f"{mw_p:.2e}".split("e")
    print(f"KS test (dr vs non {column}): statistic = {ks_res.statistic:.3f}, p-value = {ks_base} × 10^{int(ks_exp)}")
    print(f"Mann-Whitney test (dr vs non {column}): statistic = {mw_res.statistic:.3f}, p-value = {mw_base} × 10^{int(mw_exp)}")
    print()

    return fig, ax



def generate_ABC_gene_df(df_filtered_hic, ABC_threshold=0.041, hic_threshold = (-0.5), distance_threshold=50):
    # Get unique genes
    genes = df_filtered_hic['TargetGene'].unique()
    df_gene_abc = pd.DataFrame({'TargetGene': genes})
    # Per-gene calculations
    for gene in genes:
        enhancers = df_filtered_hic.loc[df_filtered_hic['TargetGene'] == gene]
        enhancers_FC = enhancers.loc[enhancers['%Change.normalized.HiC.5Kb'] <= hic_threshold ]
        enhancers_FC_ABC = enhancers_FC.loc[enhancers_FC['ABC.Score.noAux'] >= ABC_threshold]
        enhancers_distal = enhancers.loc[enhancers['DistanceToTSS.Kb'] >= distance_threshold]
        enhancers_distal_ABC = enhancers_distal.loc[enhancers_distal['ABC.Score.noAux'] >= ABC_threshold]
        total = len(enhancers)
        df_gene_abc.loc[df_gene_abc['TargetGene'] == gene, [
            'n_enhancers_distal', 'n_enhancers_proximal', 'n_enhancers_hi_HiC', 'n_enhancers_low_HiC',
            'Fraction_enhancers_distal', 'Fraction_enhancers_distal_ABC',
            'Fraction_enhancers_hic', 'Fraction_enhancers_hic_ABC',
            'n_enhancers_distal_ABC', 'n_enhancers_hic_ABC'
        ]] = [
            len(enhancers_distal),
            len(enhancers) - len(enhancers_distal),
            len(enhancers_FC),
            len(enhancers) - len(enhancers_FC),
            len(enhancers_distal) / total,
            len(enhancers_distal_ABC) / total,
            len(enhancers_FC) / total,
            len(enhancers_FC_ABC) / total,
            len(enhancers_distal_ABC),
            len(enhancers_FC_ABC)
        ]
    # CLean df
    df_gene_abc.dropna(inplace=True)
    return df_gene_abc



def subset_group_from_ABC_gene_df(df_gene_abc, df_filtered_group):
    genes = df_filtered_group['TargetGene'].unique()
    return df_gene_abc[df_gene_abc['TargetGene'].isin(genes)]
    


def calc_gene_enhancer_stats(df_gene_abc_group, group_name,
                             column,
                             column_ABC, label):
    # Compute stats
    total_genes = len(df_gene_abc_group)
    distal_count = df_gene_abc_group[df_gene_abc_group[column] >= 1].shape[0]
    strong_distal_count = df_gene_abc_group[df_gene_abc_group[column_ABC] >= 1].shape[0]
    # Print stats
    print(f"{group_name} genes:")
    print(f"Total genes: {total_genes}")
    print(f"Genes with ≥1 {label} enhancer: {distal_count}")
    print(f"Genes with ≥1 strong {label} enhancer: {strong_distal_count}")
    if total_genes > 0:
        print(f"Percent strong {label}: {strong_distal_count / total_genes:.3f}")
    else:
        print(f"Percent strong {label}: N/A (no genes)")
    print("")



def plot_fraction_cdf_and_stats(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr, fraction_col,
                                 labels=None,
                                 colors=None,
                                 threshold=0.5,
                                 save_path=None):
    # Set labels and colors
    if labels is None:
        labels = ["Cohesin-independent genes", "mid genes", "Cohesin-dependent genes"]
    if colors is None:
        colors = ['#646363', 'pink', '#F8B500']
    # Plot
    fig, ax = plt.subplots(figsize=(5,5))
    sns.ecdfplot(df_gene_abc_non[fraction_col], color=colors[0], label=labels[0], lw=4)
    sns.ecdfplot(df_gene_abc_mid[fraction_col], color=colors[1], label=labels[1], lw=4)
    sns.ecdfplot(df_gene_abc_dr[fraction_col], color=colors[2], label=labels[2], lw=4)
    # Set axes
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.set_xlabel(fraction_col)
    ax.set_ylabel("CDF")
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend()
    # Save plot
    if save_path:
        plt.savefig(save_path, format='pdf', bbox_inches='tight')
    plt.close(fig)
    # Calculate percentages below and above threshold
    def perc_below_above(df, col, thresh):
        total = len(df)
        below = (df[col] < thresh).sum()
        above = total - below
        return total, below, above
    non_total, non_below, non_above = perc_below_above(df_gene_abc_non, fraction_col, threshold)
    dr_total, dr_below, dr_above = perc_below_above(df_gene_abc_dr, fraction_col, threshold)
    # Print stats
    print(f"Non genes (total {non_total}):")
    print(f"Non genes with > {threshold*100}% {fraction_col}: n = {non_above}, fraction ={non_above / non_total:.3f}")
    print(f"DR genes (total {dr_total}):")
    print(f"DR genes with > {threshold*100}% {fraction_col}: n = {dr_above}, fraction ={dr_above / dr_total:.3f}")
    print("")
    # Odds ratio and Fisher exact test on contingency table
    contingency_table = [[dr_above, dr_below], [non_above, non_below]]
    or_res = odds_ratio(contingency_table)
    fisher_res = fisher_exact(contingency_table, alternative='two-sided')
    print(f"Odds ratio (DR above threshold vs Non above threshold): {or_res}")
    print(f"Fisher exact test result: {fisher_res}")

    return fig, ax



def plot_enhancer_count_distribution(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr, column,
                                     group_labels=None, colors=None,
                                     save_path=None):
    # Define group labels and colors
    if group_labels is None:
        group_labels = ["Cohesin-independent genes", "Cohesin-mid genes", "Cohesin-dependent genes"]
    if colors is None:
        colors = ['#646363', 'pink', '#F8B500']
    # Define the function for computing the enhnacer fractions
    def get_fractions(df):
        total = len(df)
        zero = (df[column] == 0).sum()
        one = (df[column] == 1).sum()
        two = (df[column] == 2).sum()
        three_plus = (df[column] >= 3).sum()
        return total, zero/total, one/total, two/total, three_plus/total
    # Calculate fractions for each group
    totals = []
    fractions = []
    for df in [df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr]:
        total, zero_f, one_f, two_f, three_plus_f = get_fractions(df)
        totals.append(total)
        fractions.append([zero_f, one_f, two_f, three_plus_f])
    # Arrange plot
    x = np.arange(4)
    width = 0.25
    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.bar(x - width, fractions[0], width, color=colors[0], label=group_labels[0])
    ax.bar(x, fractions[1], width, color=colors[1], label=group_labels[1])
    ax.bar(x + width, fractions[2], width, color=colors[2], label=group_labels[2])
    # Set axes
    ax.set_xlabel(f"{column}")
    ax.set_ylabel("Fraction")
    ax.set_ylim(0, 1)
    ax.set_xticks(x)
    ax.set_xticklabels(['0', '1', '2', '3+'])
    # Clean plot
    ax.legend()
    ax.spines[['right', 'top']].set_visible(False)
    # Save plot
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    # Print fraction of genes with >=1 enhancer in non and dependent groups
    print(f"Fraction with ≥1 {column} in {group_labels[0]}: {1 - fractions[0][0]:.3f}")
    print(f"Fraction with ≥1 {column} in {group_labels[2]}: {1 - fractions[2][0]:.3f}")
    # Odds ratio and Fisher test: zero vs non-zero enhancers comparing dependent vs independent
    more_zero_dependent = totals[2] - int(fractions[2][0] * totals[2])
    zero_dependent = int(fractions[2][0] * totals[2])
    more_zero_independent = totals[0] - int(fractions[0][0] * totals[0])
    zero_independent = int(fractions[0][0] * totals[0])
    contingency_table = [[more_zero_dependent, zero_dependent], [more_zero_independent, zero_independent]]
    or_res = odds_ratio(contingency_table)
    fisher_res = fisher_exact(contingency_table, alternative='two-sided')
    print(f"Odds ratio (dependent vs independent genes zero {column}): {or_res}")
    print(f"Fisher exact test result: {fisher_res}")
    print()

    return fig, ax



def subset_df_tpm(df_gene_abc_group, df_tpm, column, label=None):
    # Get genes
    genes = df_gene_abc_group['TargetGene'].unique()
    # Subset df
    df_sub = df_tpm[df_tpm['Gene Symbol'].isin(genes)][['Gene Symbol', column]].copy()
    df_sub['Group'] = label
    return df_sub


def plot_tpm_boxplot(df_tpm_non,  df_tpm_mid, df_tpm_dr, column, colors=None, save_path=None):
    # Combine & filter out zeros
    combined_tpm = pd.concat([df_tpm_non, df_tpm_mid, df_tpm_dr], ignore_index=True)
    combined_tpm = combined_tpm[combined_tpm[column] > 0]
    # Set colors
    if colors is None:
        colors = ['#646363', 'pink', '#F8B500']
    # Plot
    fig, ax = plt.subplots(figsize=(6, 5))
    ax = sns.boxplot(data=combined_tpm, x='Group', y=column, showfliers=False, palette=colors, hue = "Group", legend=False)
    sns.stripplot(data=combined_tpm, x='Group', y=column, color='black', alpha=0.4, jitter=True, ax=ax)
    # Set axes
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.set_ylabel(f'TPM ({column}, log scale)')
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Save plot
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    #plt.close(fig)


def plot_kde_by_tpm_bin(df_tpm_dr, df_filtered_DR, expr_col, hic_col,  bins=3, save_path=None):
    # Quantile-based stratification
    df_tpm_dr = df_tpm_dr.copy()
    df_tpm_dr['Expr_bin'] = pd.qcut(df_tpm_dr[expr_col], q=bins, labels=[f"Bin{i+1}" for i in range(bins)])
    # Merge data for plotting
    merged_df = df_filtered_DR.merge(
        df_tpm_dr[['Gene Symbol', 'Expr_bin', expr_col]], 
        left_on='TargetGene', right_on='Gene Symbol', 
        how='inner'
        )
    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax = sns.kdeplot(data=merged_df, x=hic_col, hue='Expr_bin', fill=True, 
                common_norm=False, alpha=0.5)
     # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    # Counts (unique Gene Symbols) and min/max per bin
    bin_stats = (
        merged_df.groupby('Expr_bin', observed=False)
        .agg(
            min_val=(expr_col, 'min'),
            max_val=(expr_col, 'max'),
            unique_genes=('Gene Symbol', 'nunique')
            )
        .reset_index()
        )   
    # Format annotations
    annotations = [
        f"{row['Expr_bin']} (n={row['unique_genes']})\nRange: {row['min_val']:.2f} - {row['max_val']:.2f}"
        for _, row in bin_stats.iterrows()
        ]
    # Add annotations outside the plot
    plt.gca().text(
        1.05, 0.5, '\n\n'.join(annotations),
        transform=plt.gca().transAxes,
        fontsize=10,
        verticalalignment='center'
        )
    # Save plot
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close(fig)


def plot_hic_vs_distance(df, enhancers, colors=None, save_path=None):
    # Set colors
    if colors is None:
        colors = {'non': '#646363', 'mid': 'pink', 'dr': '#F8B500'}
    # Remove top and bottom 1% of Hi-C change values to get rid of extremes
    lower = df['%Change.normalized.HiC.5Kb'].quantile(0.01)
    upper = df['%Change.normalized.HiC.5Kb'].quantile(0.99)
    df = df[(df['%Change.normalized.HiC.5Kb'] >= lower) & 
            (df['%Change.normalized.HiC.5Kb'] <= upper)]
    # Log-transform distance for both dfs
    df["log10_distance"] = np.log10(df['DistanceToTSS.Kb'])
    enhancers["log10_distance"] = np.log10(enhancers['DistanceToTSS.Kb'])
    # Print correlations
    r_pearson, p_pearson = stats.pearsonr(df["log10_distance"], df['%Change.normalized.HiC.5Kb'])
    print(f"Pearson's R: {r_pearson:.2f}, p = {p_pearson:.5e}")
    print(f"Number of ABC enhancers included: {len(df)}")
    # Initialize JointGrid
    g = sns.JointGrid(
        data=df, 
        x='log10_distance', 
        y='%Change.normalized.HiC.5Kb', 
        space=0,
        height=10,
        ratio=4
        )
    g.fig.set_size_inches(12, 10) 
    # Generate central full 2D KDE
    sns.kdeplot(
        data=df, 
        x='log10_distance', 
        y='%Change.normalized.HiC.5Kb', 
        fill=True, 
        cmap='rocket_r',
        ax=g.ax_joint
        )
    # Marginal KDEs by ExprGroup
    for group in df['ExprGroup'].unique():
        sub_df = df[df['ExprGroup'] == group]
        # x-density
        sns.kdeplot(
            data=sub_df,
            x='log10_distance',
            ax=g.ax_marg_x,
            color=colors.get(group, 'grey'),
            lw=2,
            label=group
            )
        # y-density
        sns.kdeplot(
            data=sub_df,
            y='%Change.normalized.HiC.5Kb',
            ax=g.ax_marg_y,
            color=colors.get(group, 'grey'),
            lw=2
            )
    # Overlay enhancer points
    g.ax_joint.scatter(
        enhancers['log10_distance'], 
        enhancers['%Change.normalized.HiC.5Kb'],
        color='black',
        s=10,
        alpha=0.7,
        label='CRUDO Enhancers'
        )
    # Legends
    g.ax_marg_x.legend(title='ExprGroup', loc='upper center', bbox_to_anchor=(1.1, 1.2))
    g.ax_joint.legend(loc='upper right')
    # Set axes
    g.ax_joint.set_xlim(-0.5, 3)
    g.ax_joint.set_ylim(-1, 0.4)
    # Save plot
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches="tight")
    plt.close()

# ------------------ Main ------------------

def main(args):

    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)

    # Load data
    df=pd.read_csv(args.predictions_noAux, sep=",")
    df_aux=pd.read_csv(args.predictions_Aux, sep=",")
    df_TPM = pd.read_csv(args.PRO_TPM, sep=None, engine="python")
    df_PROseq = pd.read_csv(args.PRO_DeSeq2, sep=None, engine="python")
    df_HK = pd.read_csv(args.housekeeping_genes, engine="python")
    df_enhancers = pd.read_csv(args.CRUDO_enhancers, engine="python")

    # Compute fold changes
    df = add_relevant_columns(df)

    # Filter
    tpm_threshold=1
    df_filtered, df_tpm_filtered, expressed_nonHK = filter_predictions(df, df_TPM, df_HK, tpm_threshold=tpm_threshold, remove_promoters=True, remove_hk=True)
    df_aux_filtered, _, _ = filter_predictions(df_aux, df_TPM, df_HK, tpm_threshold=tpm_threshold, remove_promoters=True, remove_hk=True)
    print (f'non-housekeeping genes expressed at TPM > {tpm_threshold}: {len(expressed_nonHK)}')

    # Percent retained
    df_filtered = compute_percent_retained(df_filtered, df_aux_filtered)

    # Enhancer stats
    stats_present = enhancer_stats(df_filtered)
    stats_absent = enhancer_stats(df_aux_filtered)
    print("Cohesin present:", stats_present)
    print("Cohesin absent:", stats_absent)
    print()

    #Print average change in contact for enhancer >/< thank 50Kb
    summarize_columns(df_filtered, df_filtered['DistanceToTSS.Kb'] > 50,  "%Change.normalized.HiC.5Kb",  "Distance > 50kb")
    summarize_columns(df_filtered, df_filtered['DistanceToTSS.Kb'] <= 50,  "%Change.normalized.HiC.5Kb", "Distance ≤ 50kb")

    #Print % enhancers retained for enhancer >/< thank 50Kb
    summarize_columns(df_filtered, df_filtered['DistanceToTSS.Kb'] > 50,  "PercentRetained",  "Distance > 50kb")
    summarize_columns(df_filtered, df_filtered['DistanceToTSS.Kb'] <= 50,  "PercentRetained", "Distance ≤ 50kb")

    # Compute rolling averages (Fig. 1)
    df_sorted=compute_rolling_metrics(df_filtered, sort_by='DistanceToTSS.Kb')
    fig, ax= plot_hic_contact_rolling(df_sorted, output_dir)
    save_path = os.path.join(output_dir, "ABC_dist_hic_rolling" + ".pdf")
    fig.savefig(save_path, format="pdf", bbox_inches="tight")
    fig, ax = plot_abc_score_rolling(df_sorted, output_dir)
    save_path = os.path.join(output_dir, "ABC_dist_ABC_rolling" + ".pdf")
    fig.savefig(save_path, format="pdf", bbox_inches="tight")
    plot_log2fc_hic_contact_rolling(df_sorted, output_dir, filename="ABC_dist_log2_hic_rolling")
    plot_log2fc_abc_score_rolling(df_sorted, output_dir, filename="ABC_dist_log2_ABC_rolling")
    plot_log2fc_percent_retained_rolling(df_sorted, output_dir, filename="ABC_dist_log2_retained_rolling")

    #Plot supplemental figures using ABC E-P pairs (Extended Data Fig. 1)
        # scatters
    plot_hic_contact_scatter(df_sorted, output_dir, filename="ABC_dist_hic_scatter")
    fig, axes = plot_abc_score_scatter(df_sorted, output_dir)
    save_path = os.path.join(output_dir, "ABC_dist_ABC_scatter" + ".pdf")
    fig.savefig(save_path, format="pdf", bbox_inches="tight")
    plot_log2fc_hic_contact_scatter(df_sorted, output_dir, filename="ABC_dist_log2_hic_scatter")
    plot_log2fc_abc_score_scatter(df_sorted, output_dir, filename="ABC_dist_log2_ABC_scatter")
        # histogram
    fig, ax = plot_distance_histogram(df_sorted, output_dir)
    save_path = os.path.join(output_dir, "ABC_dist_histogram" + ".pdf")
    fig.savefig(save_path, format="pdf", bbox_inches="tight")
        # CDFs
    plot_abc_enhancers_cdf(
    df_filtered, df_aux_filtered,
    gene_col='TargetGene', 
    dist_col='DistanceToTSS.Kb', 
    label_present='Cohesin present', 
    label_absent='Cohesin absent',
    color_present='#6E7CA0',
    color_absent='#8F3237',
    save_path_enhancers= os.path.join(output_dir, "ABC_n_enhancer_CDF.pdf"), 
    save_path_distance= os.path.join(output_dir, "ABC_distance_enhancers_CDF.pdf")
    )

    # Plot ABC boxplots
    plot_h3k27_hic_boxplot(df_filtered, output_dir, filename="ABC_h3k27_hic_boxplot")
    plot_retained_hic_boxplot(df_filtered, output_dir, filename="ABC_retained_hic_boxplot")
    plot_retained_dist_boxplot(df_filtered, output_dir, filename="ABC_retained_dist_boxplot")
    plot_cohesin_distance_boxplot(df_filtered, df_aux_filtered, output_dir, filename="ABC_cohesin_dist_boxplot")
    plot_CTCF_ABC_boxplot(df_filtered, output_dir, filename="ABC_CTCF_ABC_boxplot")
    plot_CTCF_hic_boxplot(df_filtered, output_dir, filename="ABC_CTCF_hic_boxplot")

    # Print # of retained vs. lost enhancers
    print(f"Enhancers retained: {len(df_filtered.loc[df_filtered['PercentRetained'] == 100])}")
    print(f"Enhancers lost: {len(df_filtered.loc[df_filtered['PercentRetained'] == 0])}")
    # Print # of enhancers with CTCF binding site clos vs without
    print(f"Enhancers with a CTCF binding site within 5Kb: {len(df_filtered.loc[df_filtered['CTCFwithin5Kb'] == True])}")
    print(f"Enhancers without a CTCF binding site within 5Kb: {len(df_filtered.loc[df_filtered['CTCFwithin5Kb'] == False])}")

    # Plot log2FC ABC score vs Hi-C
    plot_density_correlation(df_filtered, "ABC.Score.log2FC.(Aux/noAux)", 'normalized.HiC.5Kb.log2FC.(Aux/noAux)',
                         xlim=(-10,6), ylim=(-10,6), 
                         save_path=os.path.join(output_dir, "ABC_log2_ABC_HiC_density.pdf")
                         )
    ## Plots fractions of retained and lost enhancers stratified by their distance category (distal vs proximal) 
    # and proximity to CTCF (CTCF within 5Kb true or false) as bar graphs (liek Fig. 1j)
    plot_ctcf_distal_proximal_retention(df_filtered, distance_col='DistanceToTSS.Kb', ctcf_col='CTCFwithin5Kb',
                                           distal_threshold=50, 
                                           save_path=os.path.join(output_dir, "ABC_retained_CTCF_bar.pdf")
                                           )
    
    # Startify ABC-predictions based on expression category and set these thresholds
    df_filtered_hic, df_filtered_non, df_filtered_mid, df_filtered_dr = filter_proseq_and_hic(
                                                                        df_PROseq, df_TPM, df_filtered,
                                                                        TPM_threshold=5, 
                                                                        log2fc_dr_threshold=-(np.log2(1.8)),log2fc_non_threshold=-(np.log2(1.4)), 
                                                                        padj_threshold=0.05, hic_threshold=5)
    
    # Plot E-G pair comparison CDFs for the expression groupes
    # Compare distance distribution (like Extended Data Figure 15c)
    plot_cdf_comparison(df_filtered_non, df_filtered_mid, df_filtered_dr, "DistanceToTSS.Kb", 
                        xlim=None, log_scale=True, threshold=50, vline=50, 
                        save_path=os.path.join(output_dir, "ABC_ExpressionGroup_distance_CDF.pdf")
                        )
    # Compare change in 3D contact distribution (like Extended Data Figure 6c)
    plot_cdf_comparison(df_filtered_non, df_filtered_mid, df_filtered_dr, "%Change.normalized.HiC.5Kb", 
                        xlim=(-1,1), threshold=(-0.5),vline=0, 
                        save_path=os.path.join(output_dir, "ABC_ExpressionGroup_changeHiC_CDF.pdf")
                        )
    
    #G enerate per gene summary table of ABC predictions defining the threshold for distal 
    # and contact sensitive enhancers, and the ABC score that constitutes 15% effect sizes 
    # from the correlation analysis in "_06_CRUDO_Enhancer_ElementLevel_Visualization.py"   
    df_gene_abc = generate_ABC_gene_df(df_filtered_hic, ABC_threshold=0.041, hic_threshold = (-0.5), distance_threshold=50)

    # Split into expression groups
    df_gene_abc_non = subset_group_from_ABC_gene_df(df_gene_abc, df_filtered_non)
    df_gene_abc_mid = subset_group_from_ABC_gene_df(df_gene_abc, df_filtered_mid)
    df_gene_abc_dr  = subset_group_from_ABC_gene_df(df_gene_abc, df_filtered_dr)

    # Print number of distal/contact sensitive predicted enhancers for expression groups
    # Define groups for printing enhancer/gene stats
    groups = [
        (df_gene_abc_non, "Non"),
        (df_gene_abc_mid, "Mid"),
        (df_gene_abc_dr, "DR")
        ]
    # Define enhancer types: (total_col, strong_col, label)
    enhancer_types = [
            ("n_enhancers_distal", "n_enhancers_distal_ABC", "distal"),
            ("n_enhancers_hi_HiC", "n_enhancers_hic_ABC", "contact sensitive")
                ]
    # Loop over both
    for df_group, group_name in groups:
        for total_col, strong_col, enhancer_label in enhancer_types:
            calc_gene_enhancer_stats(df_group, group_name, total_col, strong_col, enhancer_label)

    # Plot per gene enhancer fraction comparison CDFs for the expression groupes,
    # specifiy the threshold for printing stats (ie 0.5 is what fraction of genes have than 50% enhancers of the specified columb )
    # Compare fraction of distal enhancers (like Extended Data Figure 15d)
    plot_fraction_cdf_and_stats(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr, "Fraction_enhancers_distal",
                                 threshold=0.5,
                                save_path=os.path.join(output_dir, "ABC_ExpressionGroup_FractionDistal_CDF.pdf")
                                )
    # Compare fraction of contact sensitive enhancers (like Fig. 6d)
    plot_fraction_cdf_and_stats(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr, "Fraction_enhancers_hic",
                                 threshold=0.5,
                                 save_path = os.path.join(output_dir, "ABC_ExpressionGroup_FractionHiC_CDF.pdf")
                                 )

    # Plot fractions of genes by number of colum-specified enhancers (0,1,2,3+), 
    # and compute odds ratio and Fisher test for zero vs non-zero enhancers between dependent and independent groups.
    # For strong distal enhancers (like Extended Data Figure 15e)
    plot_enhancer_count_distribution(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr,
                                 column='n_enhancers_distal_ABC',
                                 save_path=os.path.join(output_dir, "ABC_ExpressionGroup_nDistal_bar.pdf")
                                 )
    # For strong contact sensitive enhancers (like Figure 6e)
    plot_enhancer_count_distribution(df_gene_abc_non, df_gene_abc_mid, df_gene_abc_dr,
                                 column="n_enhancers_hic_ABC",
                                 save_path=os.path.join(output_dir, "ABC_ExpressionGroup_nHiC_bar.pdf")
                                 )
    
    # Subset TPM df based on gene expression groups to compare TPM distribution across expression groups
    df_tpm_non = subset_df_tpm(df_gene_abc_non, df_TPM, "RAD21_NT_avg", "non-DR")
    df_tpm_mid = subset_df_tpm(df_gene_abc_mid, df_TPM,  "RAD21_NT_avg","mid-DR")
    df_tpm_dr  = subset_df_tpm(df_gene_abc_dr,  df_TPM, "RAD21_NT_avg", "DR")

    # Plot Boxplot TPM comparison for expression groups
    plot_tpm_boxplot(df_tpm_non, 
                    df_tpm_mid, 
                    df_tpm_dr, "RAD21_NT_avg", 
                    save_path=os.path.join(output_dir, "ABC_ExpressionGroup_TPM_boxplot.pdf")
                    )
    
    # Plot KDE of Hi-C change stratified by quantile bins of TPM expression for strongly down-regulated genes
    plot_kde_by_tpm_bin(df_tpm_dr, df_filtered_dr, 
                    "RAD21_NT_avg", "%Change.normalized.HiC.5Kb", bins=3,
                    save_path=os.path.join(output_dir, "ABC_ExpressionGroup_DR_TPM_density.pdf")
                    )
    # Add expression group cplumn to df_filtered to prep for plotting
    df_filtered_hic['ExprGroup'] = 'non'
    df_filtered_hic.loc[df_filtered_hic['TargetGene'].isin(df_filtered_mid['TargetGene']), 'ExprGroup'] = 'mid'
    df_filtered_hic.loc[df_filtered_hic['TargetGene'].isin(df_filtered_dr['TargetGene']), 'ExprGroup'] = 'dr'

    #Filter data, calculate correlations, and plot Hi-C change vs distance with KDEs 
    # (central all, and margins stratified by group) and CRUDO enhancer points on top
    # Add expression group cplumn to df_filtered to prep for plotting
    df_filtered_hic['ExprGroup'] = 'non'
    df_filtered_hic.loc[df_filtered_hic['TargetGene'].isin(df_filtered_mid['TargetGene']), 'ExprGroup'] = 'mid'
    df_filtered_hic.loc[df_filtered_hic['TargetGene'].isin(df_filtered_dr['TargetGene']), 'ExprGroup'] = 'dr'
    #Generate the denisty plot and compute the corrleation
    plot_hic_vs_distance(df_filtered_hic, df_enhancers,
                     save_path = os.path.join(output_dir, "ABC_changeHiC_distance_density.pdf")
                     )
    
if __name__ == "__main__":
    args=parse_args()
    main(args)