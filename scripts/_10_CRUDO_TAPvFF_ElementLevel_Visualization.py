#!/usr/bin/env python3

"""
10_CRUDO_TAPvFF_ElementLevel_Visualization.py
Author: Philine Guckelberger
Date: 2025/08/07

Description:
    This script visualizes enhancer-level CRUDO effects from both FloWFISH (FF) and DC-TAP-seq (TAP) readouts.
    It generates comparative scatter and bar plots to explore correlations between the two read out methods.

Inputs:
    - A CSV file containing CRUDO FF element effects and extracted feature data (output from _04a_CRUDO_FF_datasets_merge.py)
    - A CSV file containing CRUDO DC-TAP element effects and extracted feature data (output from _09_CRUDO_TAP_analysis_pipeline.py)


Outputs:
    - Scatter plots comparing enhancer effects between FF and TAP, colored by gene
    - Bar plots comparing expression across FF and TAP for TSS elements

Usage:
    python scripts/_10_CRUDO_TAPvFF_ElementLevel_Visualization.py \
        --FF_effects resources/CRUDO_FF_SupplementaryTable2c.csv \
        --TAP_effects resources/CRUDO_TAP_ElementLevel_analysis.csv \
        --output_directory_plot path/to/output/directory/plots/

"""

# ------------------ Imports ------------------
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import stats
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Compare CRUDO element effects from FlowFISH and DC-TAP-seq.")
    parser.add_argument("--FF_effects", required=True, help="CRUDO FF element effects csv. Output from _04a_CRUDO_FF_datasets_merge.py")
    parser.add_argument("--TAP_effects", required=True, help="CRUDO TAP element effects csv. Output from _09_CRUDO_TAP_analysis_pipeline.py")
    parser.add_argument("--output_directory_plot", required=True, help="Output diretory for plots.")
    return parser.parse_args()


def grouped_scatter(df, x, y, group_col, colors, xlim=(-1,1), ylim=(-1,1),  xscale="linear", yscale="linear",
                    xticks=None, yticks=None, title=None, filename=None, output=None, yerr_col=None, xerr_col=None):
    """Scatter grouped by category with custom colors, Pearson correlation, and optional error bars."""
    fig, ax = plt.subplots(figsize=(5, 5)) 
    for group, subdf in df.groupby(group_col, observed=False):
        # Plot scatter points
        subdf.plot(
            x=x, y=y, kind='scatter', lw=1, s=125, alpha=1, legend=True,
            color=colors[group], ax=ax
        )
        # Plot error bars if yerr_col or xerr_col provided
        if yerr_col or xerr_col:
            yerr = subdf[yerr_col] if yerr_col else None
            xerr = subdf[xerr_col] if xerr_col else None
            ax.errorbar(
                subdf[x], subdf[y], 
                xerr=xerr, yerr=yerr,
                fmt='none', ecolor='gray', linewidth=0.5, zorder=-1 
            )
    # Set log scales if specified
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    # Only set limits if they are not "" (empty string)
    if xlim != "":
        ax.set_xlim(*xlim)
    if ylim != "":
        ax.set_ylim(*ylim)
    # Draw horizontal dashed line at y=0 if 0 in ylim range
    if ylim != "" and ylim[0] <= 0 <= ylim[1]:
        ax.axhline(0, linestyle='--', color='gray')
    # Draw vertical dashed line at x=0 if 0 in xlim range
    if xlim != "" and xlim[0] <= 0 <= xlim[1]:
        ax.axvline(0, linestyle='--', color='gray')
    # Clean plot
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if xticks: ax.set_xticks(xticks)
    if yticks: ax.set_yticks(yticks)
    ax.spines[['right', 'top']].set_visible(False)
    # Compute correlation
    x_vals = df[x].copy()
    y_vals = df[y].copy()
    # Add correlation line
    slope, intercept = np.polyfit(x_vals, y_vals, deg=1)
    ax.plot(x_vals, slope *  x_vals + intercept, color='black', linestyle='dashed', linewidth=0.5)
    #ste log if needed
    if xscale == "log":
        x_vals = np.log10(x_vals)
    if yscale == "log":
        y_vals = np.log10(y_vals)
    # Compute Pearson correlation on filtered & transformed values
    r, p = stats.pearsonr(x_vals, y_vals)
    if title:
        print(f"{title} R = {r:.5f}, P = {p:.3e}")
    if filename:
        path_plot = os.path.join(output, filename)
        plt.savefig(path_plot, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return (fig, ax)


def FeatureColorCorr(elements, feat_colors, group_col, output_dir, name=None):
#Compare FF vs TAP CRUDO effects
    #EnhancerEffect no Auxin
    fig, ax = grouped_scatter(elements,
            x="EnhancerEffect.noAux_FF", y='EnhancerEffect.noAux_TAP', xlim=(0,1),ylim=(0,1),
            xerr_col="ci95.EnhancerEffect.noAux_FF",yerr_col="ci95.EnhancerEffect.noAux_TAP",
            group_col=group_col, colors=feat_colors, title='Colored by {name}: EnhancerEffect noAux FF vs TAP',
            filename=f'FFvTAP_per_{name}_EnhancerEffect_noAux_scatter.pdf', output=output_dir
            )
    #EnhancerEffect plus Auxin
    fig, ax = grouped_scatter(elements,
            x="EnhancerEffect.Aux_FF", y='EnhancerEffect.Aux_TAP', xlim=(0,1),ylim=(0,1),
            xerr_col="ci95.EnhancerEffect.Aux_FF",yerr_col="ci95.EnhancerEffect.Aux_TAP",
            group_col=group_col, colors=feat_colors, title='Colored by {name}: EnhancerEffect Aux FF vs TAP',
            filename=f'FFvTAP_per_{name}_EnhancerEffect_Aux_scatter.pdf', output=output_dir
            )
    #Cohesin Dependence
    fig, ax = grouped_scatter(elements,
            x="CohesinDependence_FF", y='CohesinDependence_TAP',
            xerr_col="ci95.CohesinDependence_FF",yerr_col="ci95.CohesinDependence_TAP",
            group_col=group_col, colors=feat_colors, title='Colored by {name}: CohesinDep FF vs TAP',
            filename=f'FFvTAP_per_{name}_CohesinDep_scatter.pdf', output=output_dir
            )
    

def bar_plot(df, x, y_base, output_dir):
    df = df.loc[(df['category']=='TSS')].copy()
    # Define color and data columns
    df['FF'] = 1 - df[f'{y_base}_FF']
    df['TAP'] = 1 - df[f'{y_base}_TAP']
    color_values = {'FF': '#6D5D6E', 'TAP': '#F4EEE0'}
    # Select plotting columns
    y_cols = ['FF', 'TAP']
    # Prepare error bars
    yerr_dict = {
    'FF': df[f'ci95.{y_base}_FF'],
    'TAP': df[f'ci95.{y_base}_TAP']}
    # Plot
    fig, ax = plt.subplots(figsize=(5, 4))
    df.plot(
        x=x, y=y_cols, kind='bar', 
        color=[color_values[col] for col in y_cols], 
        legend=True, width=0.65, 
        yerr=yerr_dict, edgecolor='gray', ax=ax)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Gene expression remaining")
    ax.set_xlabel("")
    # Save plot
    path_plot = os.path.join(output_dir, f'FFvTAP_TSS_{y_base}_bar.pdf')
    plt.savefig(path_plot, format="pdf", bbox_inches="tight")
    plt.close(fig)
    
    return fig, ax




# ------------------ Main ------------------

def main(args):
    ### read in tapseq results and FF results
    df_FF =  pd.read_csv(args.FF_effects, sep=None, engine="python")
    df_TAP =  pd.read_csv(args.TAP_effects, sep=None, engine="python")

    #Merge and filter for minimum effect size
    df_combined=df_FF.merge(df_TAP, left_on=['name', 'TargetGene'],right_on=['name', 'TargetGene'], how='inner', suffixes=('_FF', '_TAP'))
    df_combined = df_combined.loc[(df_combined['EnhancerEffect.noAux_FF']>=0.05) & (df_combined['EnhancerEffect.noAux_TAP']>=0.05) ]

    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)

    #Plot CRUDO effect sizes FF vs TAP as scatters (effect no aux, effect plus aux, cohesin dep)
        #Extended Data Fig. 10a-c
    gene_colors= {'CCND1':'#85B09A', 'KITLG': '#8E9089', 'SSFA2': '#F3E600', 'FAM3C': 'purple', 'MYC': '#6E7CA0'} 
    gene_col="TargetGene"
    FeatureColorCorr(df_combined, gene_colors, gene_col, output_dir, name="Gene")

    # Plot CRUDO TSS effect comparison FF vs TAP
    #   Extended Data Fig. 10d,e
    #noAux
    bar_plot(df_combined, gene_col, "EnhancerEffect.noAux", output_dir)
    #plusAux
    bar_plot(df_combined, gene_col, "EnhancerEffect.Aux", output_dir)



if __name__ == "__main__":
    args = parse_args()
    main(args)