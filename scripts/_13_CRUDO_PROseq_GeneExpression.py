#!/usr/bin/env python3
"""
13_CRUDO_PROseq_GeneExpression.py
Author: Philine Guckelberger
Date: 2025/08/11

Description:
    This script processes PRO-seq RPKM data and converts it into TPM values.
    It computes expression percentiles for target genes and visualizes their
    expression distributions. Normalized bar plots are generated to compare
    expression with and without cohesin.


Inputs:
    - A PRO-seq RPKM txt file (in this case, GSE106886_Rao-2017-Genes.rpkm.txt)


Outputs:
    - A CSV file of TPM-normalized values
    - PDF scatter plot of percentile ranks (all genes, highlighted target genes)
    - Per-gene bar plots of normalized expression (%), with and without cohesin


Usage:
    python scripts/_13_CRUDO_PROseq_GeneExpression.py \
        --PRO_RPKM resources/GSE106886_Rao-2017-Genes.rpkm.txt \
        --output_directory path/to/output/directory/

"""


# ------------------ Imports ------------------
import argparse
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
from scipy import stats
import math
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Convert PRO-seq RPKM to TPM and visualize gene expression.")
    parser.add_argument("--PRO_RPKM", help="Path to input PRO-seq RPKM txt file")
    parser.add_argument("--output_directory", help="Directory to save outputs")
    return parser.parse_args()



def rpkm_to_tpm(rpkm):
    rpkm_sum = rpkm.sum()
    if rpkm_sum == 0:
        return rpkm
    return (rpkm / rpkm_sum) * 1e6


def calc_ci(row, cols):
    selected_cols = row[cols]
    mean = np.mean(selected_cols)
    sem = np.std(selected_cols, ddof=1) / np.sqrt(len(selected_cols))
    return 1.96 * sem


def percentile_rank(tpm_df, genes):
    # Extract gene expression values for '161219-Ctrl'
    gene_expression_values = tpm_df['RAD21_NT_avg']
    # Loop through each gene to find its percentile rank and prunt it 
    for gene in genes:
        gene_value = tpm_df.loc[tpm_df['Gene Symbol'] == gene, 'RAD21_NT_avg'].values
        if len(gene_value) > 0:
            gene_value = gene_value[0]
            # Calculate the percentile rank of the gene expression value
            percentile_rank = stats.percentileofscore(gene_expression_values, gene_value)
            print(f"{gene} is in the {percentile_rank:.2f} percentile")
        else:
            print(f"{gene} not found in the dataset")   


def plot_percentiles(tpm_df,genes, output_dir):
    # define plot
    fig = plt.figure(figsize=(10, 6))
    ax = plt.gca()
    # Plot all genes (TPM vs Percentile) with a log scale for the x-axis
    tpm_values = tpm_df['RAD21_NT_avg']
    percentiles = [stats.percentileofscore(tpm_values, tpm) for tpm in tpm_values]  
    ax.scatter(tpm_values, percentiles, color='powderblue', alpha=0.25, label='All genes')
    # Highlight the target genes
    colors = {'CCND1':'#85B09A', 'KITLG': '#8E9089', 'SSFA2': '#F3E600', 'FAM3C': 'purple', 'MYC': '#6E7CA0'}
    for gene in genes:
        gene_value = tpm_df.loc[tpm_df['Gene Symbol'] == gene, 'RAD21_NT_avg'].values
        if len(gene_value) > 0:
            # Find the percentile of this specific gene
            gene_percentile = stats.percentileofscore(tpm_values, gene_value[0])
            # Plot the gene with its assigned color from the dictionary
            ax.scatter(gene_value[0], gene_percentile, color=colors[gene], label=gene, edgecolor='grey')
    # Set axes
    ax.set_xscale('log')
    ax.set_xlim(1,) 
    # Clean plot
    plt.ylabel('Percentile Rank (%)')
    plt.xlabel('TPM')
    plt.legend()
    # save plot
    path= os.path.join(output_dir, f'CRUDO_gene_TPM_percentiles.pdf')
    fig.savefig(path, format="pdf", bbox_inches="tight")
    print(f'Saving percentiles plot to {path}')
    plt.close(fig)



def normalize(row, col, target_col):
    target_value = row[target_col]
    if target_value == 0:
        return 0
    else:
        return (row[col] / target_value) * 100
    

def map_condition(condition, grouped_cols):
    for grouped_condition, cols in grouped_cols.items():
        if condition in cols:
            return grouped_condition
    return condition



def prep_gene_expr_values_for_plotting(tpm_df):
    columns_to_norm = [
        "RAD21_NT_avg", "RAD21_T_avg", "RAD21_NT_ci", "RAD21_T_ci",
        'RAD21-mAC-1-NT-NS',  'RAD21-mAC-2-NT-S', 'RAD21-mAC-3-NT-S', 'RAD21-mAC-4-NT-S',
        'RAD21-mAC-1-T-NS', 'RAD21-mAC-2-T-S', 'RAD21-mAC-3-T-S', 'RAD21-mAC-4-T-S'
        ]
    # Normalize each column
    for col in columns_to_norm:
        tpm_df[col + '_normalized'] = tpm_df.apply(normalize, args=(col, 'RAD21_NT_avg'), axis=1)
    # Define grouped columns for mapping
    grouped_cols = {
        'Aux_normalized': ['RAD21-mAC-1-T-NS_normalized', 'RAD21-mAC-2-T-S_normalized', 'RAD21-mAC-3-T-S_normalized', 'RAD21-mAC-4-T-S_normalized'],
        'noAux_normalized': ['RAD21-mAC-1-NT-NS_normalized', 'RAD21-mAC-2-NT-S_normalized', 'RAD21-mAC-3-NT-S_normalized', 'RAD21-mAC-4-NT-S_normalized'],
        }
    # Define  columns for melting
    cols_for_melting = ['Gene Symbol'] + grouped_cols['noAux_normalized'] + grouped_cols['Aux_normalized']
    tpm_df_melted = tpm_df[cols_for_melting]

    # Melt the dataframe
    tpm_df_melted = tpm_df_melted.melt(id_vars='Gene Symbol', var_name='condition', value_name='value')
    # Apply mapping to 'condition' column, pass grouped_cols as argument
    tpm_df_melted['condition'] = tpm_df_melted['condition'].apply(lambda cond: map_condition(cond, grouped_cols))

    return tpm_df_melted



def per_gene_expr_bars(genes, tpm_df_melted, tpm_df, output_dir):
    #set plot parametes
    sns.set(rc={'figure.figsize': (2, 4)})
    sns.set_style("ticks")
    #define colors and plotting order
    color_dict = {'noAux_normalized': '#6E7CA0', 'Aux_normalized': '#8F3237'}
    hue_order = ['noAux_normalized', 'Aux_normalized']
    #Per gene
    for gene in genes:
        to_use = tpm_df_melted.loc[
            (tpm_df_melted['Gene Symbol'] == gene) &
            (tpm_df_melted['condition'].isin(hue_order))
        ]
        # Print difference in gene expression between cohesin present and bsent
        normalized_diff = 100 - tpm_df.loc[tpm_df['Gene Symbol'] == gene]['RAD21_T_avg_normalized'].values
        print(f"{gene} reduction in gene expression: {normalized_diff}")
        #Plot
        fig, ax = plt.subplots(figsize=(2,4))
        sns.barplot(
            x="Gene Symbol",
            y="value",
            hue='condition',
            data=to_use,
            capsize=0,
            errorbar=('ci', 95),
            err_kws={'linewidth': 1.5},
            palette=color_dict,
            hue_order=hue_order,
            ax=ax,
            legend=False
            )
        sns.swarmplot(
            x="Gene Symbol",
            y="value",
            hue='condition',
            data=to_use,
            color='0',
            alpha=0.7,
            size=5,
            dodge=True,
            hue_order=hue_order,
            ax=ax,
            legend=False
            )
        # Clean plot
        ax.spines[['right', 'top']].set_visible(False)
        ax.set_ylabel("Normalized Expression (%)")
        # save plot
        path= os.path.join(output_dir, f'CRUDO_{gene}_TPM_bars.pdf')
        fig.savefig(path, format="pdf", bbox_inches="tight")
        print(f'Saving {gene} plot to {path}')
        plt.close(fig)

    
# ------------------ Main ------------------

def main(args):

    # Set oupput directory
    output_dir=args.output_directory
    os.makedirs(output_dir, exist_ok=True)

    # define genes
    genes=('CCND1', 'KITLG', 'MYC', 'SSFA2', 'FAM3C')

    #Read in RPKM data with the Gene column as the index
    pro_rpkm = pd.read_csv(args.PRO_RPKM, engine="python", sep="\t", index_col=0)

    #Apply to TPM conversion to all columns (works because genes are not a column but the index) & re-set the index
    tpm_df = pro_rpkm.apply(rpkm_to_tpm, axis=0)
    tpm_df = tpm_df.reset_index()


    # Define the correct column groups based on your data
    nt_cols = ['RAD21-mAC-1-NT-NS', 'RAD21-mAC-2-NT-S', 'RAD21-mAC-3-NT-S', 'RAD21-mAC-4-NT-S']
    t_cols = ['RAD21-mAC-1-T-NS', 'RAD21-mAC-2-T-S', 'RAD21-mAC-3-T-S', 'RAD21-mAC-4-T-S']

    # Compute the average TPM across samples in each group
    tpm_df["RAD21_NT_avg"] = tpm_df[nt_cols].mean(axis=1)
    tpm_df["RAD21_T_avg"] = tpm_df[t_cols].mean(axis=1)

    # Calculate the 95% confidence interval per row using apply along axis=1
    tpm_df["RAD21_NT_ci"] = tpm_df.apply(calc_ci, args=(nt_cols,), axis=1)
    tpm_df["RAD21_T_ci"]= tpm_df.apply(calc_ci, args=(t_cols,), axis=1)

    #save TPM converted file
    path= os.path.join(output_dir, "CRUDO_Genes_Pro.tpm.txt")
    tpm_df.to_csv(path, index=False, sep="\t")

    # Filter for expressed genes
    tpm_df = tpm_df.loc[tpm_df['RAD21_NT_avg']>=1]

    # Compute and print the percentile rank of each target gene
    percentile_rank(tpm_df, genes)

    #Plot the percentile ranks of all genes highlighting the target genes
    plot_percentiles(tpm_df,genes, output_dir)
    

    # Normalize expression columns to 100 in baseline conditions and melt data for plotting
    melted_norm_tpm_df = prep_gene_expr_values_for_plotting(tpm_df)

    #Per gene plot normalized TPM expression in the presence and absence of cohesin
    per_gene_expr_bars(genes, melted_norm_tpm_df, tpm_df, output_dir)


if __name__ == "__main__":
    args=parse_args()
    main(args)