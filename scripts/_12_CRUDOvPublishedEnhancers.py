#!/usr/bin/env python3
"""
12_CRUDOvPublishedEnhancers.py
Author: Philine Guckelberger
Date: 2025/08/11

Description:
    This script processes and compares enhancer elements identified by CRUDO with 
    published enhancer datasets (Fulco et al., 2019 and an E2G CRISPRi benchmark dataset).
    It classifies enhancers based on significance and direction of gene expression change,
    computes distance to target gene TSS, and generates comparative plots of enhancer properties 
    such as distance distribution and effect sizes, as well as summary statistics.

Parameters:
    Inputs:
    - A CSV file containing CRUDO FF element effects and extracted feature data (output from _04a_CRUDO_FF_datasets_merge.py)
Note: Paths to published enhancer sets are hardcoded in the script and need to be adjusted accordingly


Outputs:
    - Scatter plots of enhancer distance vs. gene expression change for CRUDO and Fulco2019 datasets.
    - CDF plots comparing distance to TSS, enhancer effect sizes, and number of enhancers per gene 
      across CRUDO, Fulco2019, and E2G benchmark datasets.
    - Printed summary statistics comparing enhancer distributions and effect sizes across datasets.


Usage:
    python scripts/_12_CRUDOvPublishedEnhancers.py \
        --CRUDO_elements resources/CRUDO_FF_SupplementaryTable2c.csv \
        --output_directory_plot /path/to/output_dir/

"""


# ------------------ Imports ------------------
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
from scipy import stats
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Compare CRUDO enhancer effects to published enhancert sets")
    parser.add_argument("--CRUDO_elements", help="Path to input CSV file with integrated Hi-C counts, element features, and CRUDO effects from _04a_CRUDO_FF_datasets_merge.py (equivalent to SupplementaryTable2c)")
    parser.add_argument("--output_directory_plot", help="Directory to save output plots")
    return parser.parse_args()

def classify_CRUDO_significant_elements(df):
    #Add Delta aux column
    df['delta.noAux'] = abs(df['EnhancerEffect.noAux.Rep1'] - df['EnhancerEffect.noAux.Rep2'])
    #Get absolute enhancer effetcts to classify both signficantly up-and down-regulated elements
    for i in range(1, 5):
        df[f'abs_EnhancerEffect.noAux.Rep{i}'] = 1 - df[f'EnhancerEffect.noAux.Rep{i}']
    #Define signficance conditions
    median_h3k27 = np.median(df['H3K27ac.RPM.values.noAux'].dropna())
    sig_conditions = (
          (df['adj.pval.EnhancerEffect.noAux.Rep1'] < 0.05) & (df['adj.pval.EnhancerEffect.noAux.Rep2'] < 0.05) &
         (df['abs_EnhancerEffect.noAux.Rep1'] >= 0.05) & (df['abs_EnhancerEffect.noAux.Rep2'] >= 0.05) &
        (df['H3K27ac.RPM.values.noAux'] > median_h3k27) &
        (df['delta.noAux'] <= 0.05) &
        (df['n'] >= 10) &
        (
        (df['TargetGene'] != 'KITLG') |
        (
            (df['TargetGene'] == 'KITLG') &
            (df['adj.pval.EnhancerEffect.noAux.Rep3'] < 0.05) &
            (df['adj.pval.EnhancerEffect.noAux.Rep4'] < 0.05) &
            (df['abs_EnhancerEffect.noAux.Rep3'] >= 0.05) &
            (df['abs_EnhancerEffect.noAux.Rep4'] >= 0.05)
        )
         )
    )
    # Define up-regulated, down-regulated, and non-significant based on Avg effect no auxin values
    df['SigCategory'] = 'non-significant'  # Default to non-significant
    df.loc[sig_conditions & (df['EnhancerEffect.noAux'] < (-0.05)), 'SigCategory'] = 'up-regulated'
    df.loc[sig_conditions & (df['EnhancerEffect.noAux'] > 0.05), 'SigCategory'] = 'down-regulated'
    # Add column for plotting average gene expression remaining
    df['Fraction change in gene expr'] = -(df['EnhancerEffect.noAux'])

    return df


def EnhanceerScatter_SigColored(for_scatter_nonsig, for_scatter_sig, output_dir, file_name=None):
    # Define colors
    colors= {'up-regulated':'#C6BF69', 'down-regulated': '#778F78', 'non-significant': '#AEA7A0'} 
    # Define plot
    fig = plt.figure(figsize=(5, 5))
    ax = plt.gca()
    # Plot non-signficiant elements in the background
    x = for_scatter_nonsig['DistanceToTSS.Kb']
    y = for_scatter_nonsig['Fraction change in gene expr']
    ax.scatter(x.values, y.values, facecolors='none', edgecolors=['#AEA7A0'], alpha=1, s=50, lw=2)
    # Plot signficiant elements on top
    x2 = for_scatter_sig['DistanceToTSS.Kb']
    y2 = for_scatter_sig['Fraction change in gene expr']
    ax.scatter(for_scatter_sig['DistanceToTSS.Kb'], for_scatter_sig['Fraction change in gene expr'], facecolors='none', 
           edgecolors=for_scatter_sig['SigCategory'].apply(lambda x: colors[x]), alpha= 1, s=50, lw=2)
    #Clean plot
    ax.set_xscale('log') 
    ax.set_xlim(0.9,10000) 
    ax.set_ylim(-1,0.5)     
    ax.axvline(x=50, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1)
    ax.axhline(y=0, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1) 
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))    
    ax.spines[['top','right']].set_visible(False)
    ax.set_xlabel('DistanceToTSS.Kb')
    ax.set_ylabel("Reduction in gene expression")
    #Set labels
    fake_handles = [mpatches.Patch(color=item) for item in colors.values()]
    label = colors.keys()
    ax.legend(fake_handles, label, loc='best', frameon=False)
    #save plot
    if file_name:
        fig.savefig(os.path.join(output_dir, file_name))
    plt.close(fig)
            
    return fig, ax
    
def PrepFulco2019Data():
    #Read in Fulco2019 CRISPRi FF data (Supp Table 3a)
    fulco2019_3a= pd.read_csv('resources/Fulco2019_S3a.txt', 
                sep=None, engine="python")
    #Use Fulco 2019 supplemetary Table 6a to map TSS positions to compute distance
    fulco2019_6a = pd.read_csv('resources/Fulco2019_S6a.txt', 
                sep=None, engine="python")
    #extract TSS position to com[pute distance
    fulco2019_genes = fulco2019_3a['Gene'].unique()
    # Filter and get only unique gene-to-TSS mappings
    fulco2019_tss = fulco2019_6a[fulco2019_6a['Gene'].isin(fulco2019_genes)][['Gene', 'Gene TSS']].drop_duplicates()
    # Convert to dictionary
    tss_dict = dict(zip(fulco2019_tss['Gene'], fulco2019_tss['Gene TSS']))
    # Add TargetGeneTSS column by mapping from dictionary
    fulco2019_3a['TargetGeneTSS'] = fulco2019_3a['Gene'].map(tss_dict)
    # Compute distance
    fulco2019_3a['enhancer_midpoint'] = ((fulco2019_3a['start'] + fulco2019_3a['end']) // 2)
    fulco2019_3a['DistanceToTSS.Kb'] = ((fulco2019_3a['enhancer_midpoint'] - fulco2019_3a['TargetGeneTSS']).abs()) / 1000
    # Classify elements
    fulco2019_3a['name'] = fulco2019_3a['Element name']+'|'+fulco2019_3a['Gene']
    fulco2019_elements = fulco2019_3a['name'].unique()
    for e in fulco2019_elements:
        if fulco2019_3a.loc[fulco2019_3a['name']== e, 'Significant'].values[0] == True:
            if ((fulco2019_3a.loc[fulco2019_3a['name'] == e, 'Fraction change in gene expr'].values[0] < -0.05) and fulco2019_3a.loc[fulco2019_3a['name'] == e, 'Valid E-G connection'].values[0] == True):
                fulco2019_3a.loc[fulco2019_3a['name']== e, 'SigCategory'] = 'down-regulated'
            elif ((fulco2019_3a.loc[fulco2019_3a['name'] == e, 'Fraction change in gene expr'].values[0] > 0.05) and fulco2019_3a.loc[fulco2019_3a['name'] == e, 'Valid E-G connection'].values[0] == True):
                fulco2019_3a.loc[fulco2019_3a['name']== e, 'SigCategory'] = 'up-regulated'
            else: 
                fulco2019_3a.loc[fulco2019_3a['name']== e, 'SigCategory'] = 'non-significant'
        else:
            fulco2019_3a.loc[fulco2019_3a['name']== e, 'SigCategory'] = 'non-significant'
    # Rename effect size column to match other data sets
    fulco2019_3a.rename(columns={"Gene":"TargetGene"}, inplace=True)
    # Stratify by significant and non-signficant elements to plot element distribution
    fulco2019_sig = fulco2019_3a.loc[(fulco2019_3a['Significant']==True) & 
                                  (fulco2019_3a['Invalidating reason']!='Enhancer overlaps bounds of gene (i.e., is promoter or intronic enhancer)')
                              ]
    fulco2019_nonsig = fulco2019_3a.loc[(fulco2019_3a['Significant']==False) &
                                     (fulco2019_3a['Invalidating reason']!='Enhancer overlaps bounds of gene (i.e., is promoter or intronic enhancer)')]
    
    return fulco2019_6a, fulco2019_3a, fulco2019_sig, fulco2019_nonsig



def CDF_2_datasets(fulco19_enhancers, CRUDO_enhancers, column, output_dir):
    # Define plot
    fig, ax = plt.subplots(figsize=(5,5))
    # Plot both data set CDFS
    sns.ecdfplot(fulco19_enhancers[column], color='#5F8D4E', label="Fulco2019 enhancers", lw=4)
    sns.ecdfplot(CRUDO_enhancers[column], color='#F7A4A4', label="CRUDO enhancers", lw=4)
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend(prop={'size': 8})
    ax.set_xlabel(column)
    ax.set_ylabel("CDF")
    #Set specific rules for distance
    if column =="DistanceToTSS.Kb":
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.axvline(x=50, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1) 
    # Print u-test statistics
    stat, p = stats.mannwhitneyu(fulco19_enhancers[column], CRUDO_enhancers[column])
    print(f"{column} U-test = {stat:.1f}, P = {p:.3e}")
    # Save figure
    filename = f"CRUDO_v_Fulco2019_CDF_{column}.pdf"
    fig.savefig(os.path.join(output_dir, filename))
    plt.close(fig)

    return fig,ax




def PrepE2GBenchmarkData():
    # Read in data
    e2g_benchmark= pd.read_csv('path/to/Combined CRISPR data - EPCrisprBenchmark_CombinedData_GRCh38.csv', sep=None, engine="python")
    # Compute distance
    e2g_benchmark['enhancer_midpoint'] = ((e2g_benchmark['chromStart']+e2g_benchmark['chromEnd'])//2)
    e2g_benchmark['TargetGeneTSS'] = ((e2g_benchmark['startTSS']+e2g_benchmark['endTSS'])//2)
    e2g_benchmark['DistanceToTSS.Kb'] = ((e2g_benchmark['enhancer_midpoint'] - e2g_benchmark['TargetGeneTSS']).abs()) / 1000
    e2g_benchmark['Reference'].unique()
    # CLassify elements
    e2g_enhancers= e2g_benchmark.loc[(e2g_benchmark['Significant']==True) &(e2g_benchmark['Regulated']==True) 
                           &(e2g_benchmark['EffectSize']<=-0.05) 
                          & (e2g_benchmark['Reference']!='Fulco et al., 2019') # Exclude Fulco2019 because we have this as a seperate data set
                          & (e2g_benchmark['Reference']!='Fulco et al., 2016') # Exclude Fulco 2016 because this covers MYC which we expect to also be cohein-sensitive in K562s
                           ].copy() 
    # Rename effect size column to match other data sets
    e2g_enhancers.rename(columns={"EffectSize": "Fraction change in gene expr", "measuredGeneSymbol":"TargetGene"}, inplace=True)

    return e2g_enhancers


def CDF_3_datasets(e2g_enhancers, fulco19_enhancers, CRUDO_enhancers, column, output_dir):
    # Define plot
    fig, ax = plt.subplots(figsize=(5,5))
    # Plot all three CDFS
    sns.ecdfplot(e2g_enhancers[column], color='gold', label="E2G enhancers", lw=4)
    sns.ecdfplot(fulco19_enhancers[column], color='darkorange', label="Fulco2019 enhancers", lw=4)
    sns.ecdfplot(CRUDO_enhancers[column], color='m',label="CRUDO enhancers", lw=4)
    # clean plot
    ax.spines[['right', 'top']].set_visible(False)
    ax.legend(prop={'size': 8})
    ax.set_xlabel(column)
    ax.set_ylabel("CDF")
    #Set specific rules for distance
    if column =="DistanceToTSS.Kb":
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        ax.axvline(x=50, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1) 
    # Save figure
    filename = f"CRUDO_v_Fulco2019_v_E2G_CDF_{column}.pdf"
    fig.savefig(os.path.join(output_dir, filename))

    plt.close(fig)

    return fig,ax


# ------------------ Main ------------------

def main(args):

    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)

    # Main
    df_CRUDO= pd.read_csv(args.CRUDO_elements, sep=None, engine="python")

    #Classify elements as significantly up-regulated, significantly down-regulated, or non-significant
    df_CRUDO_classified= classify_CRUDO_significant_elements(df_CRUDO)
    #Prepare data for plotting
    exclude_cats = ['TSS', 'negative_control']
    CRUDO_base_mask = ~df_CRUDO_classified['category'].isin(exclude_cats)
    for_scatter_sig = df_CRUDO_classified.loc[(df_CRUDO_classified['SigCategory'] != 'non-significant') & CRUDO_base_mask]
    for_scatter_nonsig = df_CRUDO_classified.loc[(df_CRUDO_classified['SigCategory'] == 'non-significant') & CRUDO_base_mask]

    #Plot distance distribution of classified elements
    EnhanceerScatter_SigColored(for_scatter_nonsig, for_scatter_sig, output_dir, "CRUDO_elements_scatter.pdf")

    #Get Fulc 2019 data and plot distribution
    fulco2019_compilation, fulco2019, fulco2019_sig, fulco2019_nonsig = PrepFulco2019Data()
    EnhanceerScatter_SigColored(fulco2019_nonsig, fulco2019_sig, output_dir, "Fulco2019_elements_scatter.pdf")


    #Plot CRUDO vs. Fulco 2010 enhancers comparisons CDFs

    #Distance distribution of enhancers
    CRUDO_enhancers=for_scatter_sig.loc[for_scatter_sig['SigCategory']=='down-regulated']
    fulco19_enhancers=fulco2019_sig.loc[fulco2019_sig['SigCategory']=='down-regulated']
    CDF_2_datasets(fulco19_enhancers, CRUDO_enhancers, "DistanceToTSS.Kb",output_dir)
    #distal enhancer effect sizes
    CRUDO_distal= CRUDO_enhancers.loc[CRUDO_enhancers['DistanceToTSS.Kb']>=50]
    fulco19_distal= fulco19_enhancers.loc[fulco19_enhancers['DistanceToTSS.Kb']>=50]
    CDF_2_datasets(fulco19_distal,CRUDO_distal,  "Fraction change in gene expr", output_dir)
    #Number of enhancers per gene
    fulco2019_genes_df = (
        fulco19_enhancers.groupby('TargetGene')
        .size()
        .reset_index(name='n_enhancers')
            )
    CRUDO_genes_df = (
        CRUDO_enhancers.groupby('TargetGene')
        .size()
        .reset_index(name='n_enhancers')
            )
    CDF_2_datasets(fulco2019_genes_df, CRUDO_genes_df,  "n_enhancers", output_dir)


    #Compare to larger E2G benchmarking CRISPRi data set

    #Distance distribution of enhancers
    e2g_enhancers = PrepE2GBenchmarkData()
    CDF_3_datasets(e2g_enhancers, fulco19_enhancers, CRUDO_enhancers, "DistanceToTSS.Kb", output_dir)
    #distal enhancer effect sizes
    e2g_distal= e2g_enhancers.loc[e2g_enhancers['DistanceToTSS.Kb']>=50]
    CDF_3_datasets(e2g_distal, fulco19_distal,CRUDO_distal,  "Fraction change in gene expr", output_dir)

    #Print dataset stats
    datasets = [
        ("E2GBenchmark", e2g_enhancers),
        ("Fulco2019", fulco19_enhancers),
        ("CRUDO", CRUDO_enhancers)
        ]
    for name, dataset in datasets:
        genes = dataset['TargetGene'].unique()
        epg = len(dataset) / len(genes)
        distal = dataset.loc[dataset['DistanceToTSS.Kb'] >= 50]
        percent_dist = len(distal) / len(dataset)
        distal_effect = distal.loc[distal['Fraction change in gene expr'] <= (-0.2)]
        print(f'{name} # of target genes: {len(genes)}')
        print(f'{name} identified enhancers: {len(dataset)}')
        print(f'{name} enhancers per gene: {epg}')
        print(f'{name} identified distal enhancers: {len(distal)}')
        print(f'{name} percent distal enhancers: {percent_dist}')
        print(f'{name} identified distal enhancers > 20% effect: {len(distal_effect)}')
        print("")


    #Print median enhancer effect sizes from data set compilation

    #Fulco2019 compilation:
    fulco2019_compilation_median_effect= np.median(
        fulco2019_compilation.loc[(fulco2019_compilation['Significant']==True)]['Fraction change in gene expr'])

    #E2G benchmark compilation:
    e2g_compilation_median_effect = np.median(e2g_enhancers['Fraction change in gene expr'])

    print(f'Compiled dataset median enhancer effects: Fulco2019 {fulco2019_compilation_median_effect} & E2G benchmark {e2g_compilation_median_effect} ')


if __name__ == "__main__":
    args=parse_args()
    main(args)