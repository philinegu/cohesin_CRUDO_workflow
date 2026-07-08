#!/usr/bin/env python3
"""
08_v2_CRUDO_FF_GuideLevel_Visualization.py
Author: Philine Guckelberger
Date: 2025/08/07 (v2)


Description:
    Script to visualize enhancer and guide level effects of CRUDO-FF data
    comparing effects in Auxin treated vs untreated conditions.
    It aggregates per-guide data across replicates, computes guide effect differences,
    and plots enhancer-level changes in gene expression.

    v2 adds a cohesin-dependence bar plot with per-gRNA dots overlaid:
      * PlotCohesinDepWithGuides -- element-level CohesinDependence bars colored by
        PowerCategory, with per-gRNA CohesinDependence.perGuide dots overlaid.
      * Bars ordered by signed distance to TSS (upstream -> TSS -> downstream).
      * Per-guide CohesinDependence.perGuide is computed inside
        create_enhancer_to_guide_dict as Delta_aux / (1 - mean(mleAvg_noAux)) per
        element, mirroring the element-level formula in _03 at gRNA granularity.


Inputs:
    - A CSV file with CRUDO enhancer regions (output from _05b_CRUDO_RAD21_analysis.py).
    - Per-guide FF outputs for each replicate (individual files per gene and replicate)

Outputs (per gene):
    - CRUDO_<gene>_Enhancer_GuideLevel_Effects.pdf      (original guide-effect plot)
    - CRUDO_<gene>_Enhancer_CohesinDep_GuideOverlay.pdf (new bars + dot overlay)


Usage:
    python scripts/_08_v2_CRUDO_FF_GuideLevel_Visualization.py  \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --processed_FF_directory resources/byExperimentRep/ \
        --output_directory_plot path/to/output/directory/plots/
"""

# ------------------ Imports ------------------

import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
plt.rcParams['pdf.fonttype'] = 42
# ------------------ HELPER FUNCTIONS ------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze CRUDO FF data and save element level data.")
    parser.add_argument("--enhancers", required=True, help="CRUDO enhancer csv. ie output from _05b_CRUDO_RAD21_analyzsis.py")
    parser.add_argument("--processed_FF_directory", required=True, help="Directory that leads to the PerExperimentReps, individual gRNA level data per gene, replicate, and condition")
    parser.add_argument("--output_directory_plot", required=True, help="Directory to save output CSV")
    return parser.parse_args()


def _signed_dist_kb(name, gene, tss_map):
    """Signed distance (Kb) enhancer midpoint -> gene TSS midpoint.
    Negative = upstream, positive = downstream. Takes a {gene: TSS_midpoint_bp} dict
    so the caller controls the coordinate set."""
    if not isinstance(name, str) or ':' not in name:
        return float('nan')
    tss_mid = tss_map.get(gene)
    if tss_mid is None:
        return float('nan')
    try:
        _, se = name.split(':')
        s, e = se.split('-')
        return ((int(s) + int(e)) / 2 - tss_mid) / 1000.0
    except Exception:
        return float('nan')


def load_df_guides_noAux(processed_FF_directory):
    file_path=processed_FF_directory
    #if all genes have the same # of replicates this can be change to pilot_genes and the KITLG part can be removed from the script
    pilot_genes_no_kit= ('CCND1', 'FAM3C', 'MYC', 'SSFA2')
    df_guides_noAux= pd.DataFrame()
    ident_columns= ['chr', 'start', 'end', 'name', 'score', 'strand', 'GuideSequence',
           'target', 'OffTargetScore', 'OligoID']

    for gene in pilot_genes_no_kit:
        replicates = [1, 2]
        dfs = []
        for rep in replicates:
            path = os.path.join(file_path, f'{gene}-0-Rep1-{rep}.scaled.txt')
            dfs.append(pd.read_csv(path, sep=None, engine="python"))
        temp_df = dfs[0].merge(dfs[1], on=ident_columns, how='inner', suffixes=('_rep1', '_rep2'))
        temp_df['target_gene']= gene
        df_guides_noAux=pd.concat([df_guides_noAux, temp_df], ignore_index=True)

    df_guides_noAux['mleAvg']=(df_guides_noAux['mleAvg_rep1']+df_guides_noAux['mleAvg_rep2'])/2
    df_guides_noAux['sumcells']= df_guides_noAux['sumcells_rep1']+df_guides_noAux['sumcells_rep2']

    #special handling for kitlg
    only_kit= ('KITLG',)
    for gene in only_kit:
        replicates = [1, 2, 3, 4]
        dfs = []
        for rep in replicates:
            path = os.path.join(file_path, f'{gene}-0-Rep1-{rep}.scaled.txt')
            dfs.append(pd.read_csv(path, sep=None, engine="python"))
        # Merge first two replicates with suffixes for rep1 and rep2
        temp_df = dfs[0].merge(dfs[1], on=ident_columns, how='inner', suffixes=('_rep1', '_rep2'))
        # Merge the third replicate (no suffix yet because now column anmes are unique after passing the previous suffixes)
        temp_df = temp_df.merge(dfs[2], on=ident_columns, how='inner')
        # Merge the fourth replicate with suffixes for rep3 and rep4 columns to finish
        temp_df = temp_df.merge(dfs[3], on=ident_columns, how='inner', suffixes=('_rep3', '_rep4'))

        temp_df['target_gene']= 'KITLG'

        temp_df['mleAvg']=(temp_df['mleAvg_rep1']+temp_df['mleAvg_rep2']+ temp_df['mleAvg_rep3']+temp_df['mleAvg_rep4'])/4
        temp_df['sumcells']= temp_df['sumcells_rep1']+temp_df['sumcells_rep2']+temp_df['sumcells_rep3']+temp_df['sumcells_rep4']

        df_guides_noAux=pd.concat([df_guides_noAux, temp_df], ignore_index=True)

    df_guides_noAux['start']=df_guides_noAux['start'].fillna(0)
    df_guides_noAux['end']=df_guides_noAux['end'].fillna(0)
    df_guides_noAux['chr']=df_guides_noAux['chr'].fillna('chr0')
    df_guides_noAux = df_guides_noAux.drop_duplicates()

    print(df_guides_noAux.columns)
    return df_guides_noAux

def load_df_guides_plusAux(processed_FF_directory):
    file_path=processed_FF_directory
    #if all genes have the same # of replicates this can be change to pilot_genes and the KITLG part can be removed from the script
    pilot_genes_no_kit= ('CCND1', 'FAM3C', 'MYC', 'SSFA2')
    df_guides_plusAux= pd.DataFrame()
    ident_columns= ['chr', 'start', 'end', 'name', 'score', 'strand', 'GuideSequence',
           'target', 'OffTargetScore', 'OligoID']

    for gene in pilot_genes_no_kit:
        replicates = [1, 2]
        dfs = []
        for rep in replicates:
            path = os.path.join(file_path, f'{gene}-6-Rep1-{rep}.scaled.txt')
            dfs.append(pd.read_csv(path, sep=None, engine="python"))
        temp_df = dfs[0].merge(dfs[1], on=ident_columns, how='inner', suffixes=('_rep1', '_rep2'))
        temp_df['target_gene']= gene
        df_guides_plusAux=pd.concat([df_guides_plusAux, temp_df], ignore_index=True)

    df_guides_plusAux['mleAvg']=(df_guides_plusAux['mleAvg_rep1']+df_guides_plusAux['mleAvg_rep2'])/2
    df_guides_plusAux['sumcells']= df_guides_plusAux['sumcells_rep1']+df_guides_plusAux['sumcells_rep2']

    #special handling for kitlg
    only_kit= ('KITLG',)
    for gene in only_kit:
        replicates = [1, 2, 3, 4]
        dfs = []
        for rep in replicates:
            path = os.path.join(file_path, f'{gene}-6-Rep1-{rep}.scaled.txt')
            dfs.append(pd.read_csv(path, sep=None, engine="python"))
        # Merge first two replicates with suffixes for rep1 and rep2
        temp_df = dfs[0].merge(dfs[1], on=ident_columns, how='inner', suffixes=('_rep1', '_rep2'))
        # Merge the third replicate (no suffix yet because now column anmes are unique after passing the previous suffixes)
        temp_df = temp_df.merge(dfs[2], on=ident_columns, how='inner')
        # Merge the fourth replicate with suffixes for rep3 and rep4 columns to finish
        temp_df = temp_df.merge(dfs[3], on=ident_columns, how='inner', suffixes=('_rep3', '_rep4'))

        temp_df['target_gene']= 'KITLG'
        temp_df['mleAvg']=(temp_df['mleAvg_rep1']+temp_df['mleAvg_rep2']+ temp_df['mleAvg_rep3']+temp_df['mleAvg_rep4'])/4
        temp_df['sumcells']= temp_df['sumcells_rep1']+temp_df['sumcells_rep2']+temp_df['sumcells_rep3']+temp_df['sumcells_rep4']

        df_guides_plusAux=pd.concat([df_guides_plusAux, temp_df], ignore_index=True)

    df_guides_plusAux['start']=df_guides_plusAux['start'].fillna(0)
    df_guides_plusAux['end']=df_guides_plusAux['end'].fillna(0)
    df_guides_plusAux['chr']=df_guides_plusAux['chr'].fillna('chr0')
    df_guides_plusAux = df_guides_plusAux.drop_duplicates()

    print(df_guides_plusAux.columns)
    return df_guides_plusAux


def create_enhancer_to_guide_dict(df_TargetList, df_guides_noAux, df_guides_plusAux):
    elements = df_TargetList['name'].unique()
    element_to_guide_dict = {}
    ident_columns= ['chr', 'start', 'end', 'name', 'score', 'strand', 'GuideSequence',
           'target', 'OffTargetScore', 'OligoID', 'target_gene']

    for element in elements:
        start = df_TargetList.loc[df_TargetList['name']==element,'start'].values[0]
        end = df_TargetList.loc[df_TargetList['name']==element,'end'].values[0]
        target_gene = df_TargetList.loc[df_TargetList['name']==element,'TargetGene'].values[0]
        df_noAux= df_guides_noAux.loc[(df_guides_noAux['start']>=start) & (df_guides_noAux['end'] <= end)
                       & (df_guides_noAux['target_gene']==target_gene)]

        df_plusAux= df_guides_plusAux.loc[(df_guides_plusAux['start']>=start) & (df_guides_plusAux['end'] <= end)
                       & (df_guides_plusAux['target_gene']==target_gene)]

        df_merge = df_noAux.merge(df_plusAux, how='inner', on=ident_columns, suffixes=('_noAux', '_plusAux'))
        df_merge['Delta_aux']= (1-(df_merge['mleAvg_noAux'])) - (1-(df_merge['mleAvg_plusAux']))
        df_merge['Delta_aux_rep1']= (1-(df_merge['mleAvg_rep1_noAux'])) - (1-(df_merge['mleAvg_rep1_plusAux']))
        df_merge['Delta_aux_rep2']= (1-(df_merge['mleAvg_rep2_noAux'])) - (1-(df_merge['mleAvg_rep2_plusAux']))
        df_merge['Delta_aux_rep3']= (1-(df_merge['mleAvg_rep3_noAux'])) - (1-(df_merge['mleAvg_rep3_plusAux']))
        df_merge['Delta_aux_rep4']= (1-(df_merge['mleAvg_rep4_noAux'])) - (1-(df_merge['mleAvg_rep4_plusAux']))

        # Compute guide effects individually
        for i in range(1, 5):
            df_merge[f'GuideEffect.noAux.Rep{i}'] = 1 - df_merge[f'mleAvg_rep{i}_noAux']
            df_merge[f'GuideEffect.Aux.Rep{i}'] = 1 - df_merge[f'mleAvg_rep{i}_plusAux']

        df_merge['GuideEffect.noAux'] = 1 - df_merge['mleAvg_noAux']
        df_merge['GuideEffect.Aux'] = 1 - df_merge['mleAvg_plusAux']

        # Per-guide CohesinDependence for the dot overlay in the bar plot below.
        # Mirrors the element-level formula in _03 at per-gRNA granularity.
        avg_guide = 1 - np.mean(df_merge['mleAvg_noAux'])
        if avg_guide != 0:
            df_merge['CohesinDependence.perGuide'] = df_merge['Delta_aux'] / avg_guide

        # Clip extreme values to stay within bounds
        effect_columns = [col for col in df_merge.columns if 'GuideEffect' in col]
        df_merge[effect_columns] = df_merge[effect_columns].clip(-1, 1)

        if element not in element_to_guide_dict.keys():
            element_to_guide_dict[element] = df_merge.copy()

    print(element_to_guide_dict[element].columns)
    return element_to_guide_dict


def PlotGuideEffects(enhancers,element_to_guide_dict, gene):
    gene_enhancers= enhancers.loc[(enhancers['TargetGene']==gene)]

    # Concatenate all guides for enhancers of the gene
    list_enhancers= gene_enhancers['name'].unique()
    dfs_to_append = []
    for enhancer in list_enhancers:
        temp_df=element_to_guide_dict[enhancer].copy()
        temp_df['enhancer']=enhancer
        dfs_to_append.append(temp_df)
        gene_enhancer_guides = pd.concat(dfs_to_append, ignore_index=True)

    # Print adjusted p-values for enhancer significance
    print(f'Adj.CohesinDependence_{gene}:')
    print(gene_enhancers['adj.pval.CohesinDependence'].map(lambda x: f"{x:.4f}"))


    # Prepare date for plotting
    data = pd.DataFrame({
    'name': gene_enhancer_guides['enhancer'],
    'noAux': gene_enhancer_guides['GuideEffect.noAux'],
    'plusAux': gene_enhancer_guides['GuideEffect.Aux']
    })
    # Melt the DataFrame to combine 'rep1' and 'rep2'
    data_melted = data.melt(id_vars=['name'], var_name='Condition', value_name='Value')

    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 6))
    # Assign x-positions to enhancers
    unique_names = data['name'].unique()
    positions = {name: i for i, name in enumerate(unique_names)}
    # Plot paired data points and their means with error bars
    for name in unique_names:
        name_data = data[data['name'] == name]
        x_position = positions[name]
        rep1_values = name_data['noAux'].values
        rep2_values = name_data['plusAux'].values
        x_rep1 = [x_position - 0.2] * len(rep1_values)
        x_rep2 = [x_position + 0.2] * len(rep2_values)
        # Scatter and connecting lines
        ax.plot(x_rep1, rep1_values, 'o', markersize=3, alpha=0.25, color="0")
        ax.plot(x_rep2, rep2_values, 'o', markersize=3, alpha=0.25, color="0")
        for x1, x2, y1, y2 in zip(x_rep1, x_rep2, rep1_values, rep2_values):
            ax.plot([x1, x2], [y1, y2], 'k-', linewidth=0.25, alpha=0.25)
        # Plot means and confidence intervals
        mean_rep1_values = np.mean(name_data['noAux'].values)
        mean_rep2_values = np.mean(name_data['plusAux'].values)
        x_mean_rep1 = [x_position-0.2]
        x_mean_rep2 = [x_position+0.2]
        # Errors
        stdev_noAux = np.std(name_data['noAux'],ddof=1)
        stdev_plusAux = np.std(name_data['plusAux'],ddof=1)
        sqrt_n = math.sqrt(len(element_to_guide_dict[enhancer]))
        mean_rep1_ci = 1.96*(stdev_noAux/sqrt_n)
        mean_rep2_ci= 1.96*(stdev_plusAux/sqrt_n)
        plt.errorbar(x_mean_rep1, mean_rep1_values, yerr=mean_rep1_ci, fmt='_', label='No Auxin', color='#6E7CA0', alpha=1, markersize=25)
        plt.errorbar(x_mean_rep2, mean_rep2_values, yerr=mean_rep2_ci, fmt='_', label='Plus Auxin Rep', color='#8F3237', alpha=1, markersize=25)

    # Clean plot
    ax.set_xticks([positions[name] for name in unique_names])
    x_labels = gene_enhancers['DistanceToTSS.Kb'].unique()
    ax.set_xticklabels(x_labels, rotation=0, ha="right")
    ax.set_ylabel("Reduction in gene expression")
    ax.spines[['right', 'top']].set_visible(False)
    ax.axhline(0, linestyle='--', color='gray')

    plt.close()

    return fig,ax


def PlotCohesinDepWithGuides(enhancers, element_to_guide_dict, gene,
                             cohesin_dep_colors, gene_tss_midpoint,
                             dot_color='0', dot_alpha=0.25, dot_size=2,
                             ylim=(-0.55, 1)):
    """Element-level CohesinDependence bars colored by PowerCategory, with per-gRNA
    CohesinDependence.perGuide dots overlaid. Bars ordered by signed distance to TSS
    (upstream -> TSS -> downstream)."""
    gene_enh = enhancers.loc[enhancers['TargetGene'] == gene].copy()
    gene_enh['signed_distance_kb'] = [
        _signed_dist_kb(n, gene, gene_tss_midpoint) for n in gene_enh['name']
    ]
    gene_enh = gene_enh.sort_values('signed_distance_kb', kind='stable',
                                    na_position='last').reset_index(drop=True)
    unique_names = list(gene_enh['name'])

    fig, ax = plt.subplots(figsize=(max(6, 1.2 * len(unique_names)), 5))
    positions = {n: i for i, n in enumerate(unique_names)}

    bar_x = [positions[n] for n in unique_names]
    bar_y = gene_enh['CohesinDependence'].values
    bar_err = (gene_enh['ci95.CohesinDependence'].values
               if 'ci95.CohesinDependence' in gene_enh.columns else None)
    if 'PowerCategory' in gene_enh.columns:
        bar_colors = gene_enh['PowerCategory'].map(cohesin_dep_colors).fillna('lightgray').values
    else:
        bar_colors = ['lightgray'] * len(bar_x)
    ax.bar(bar_x, bar_y, width=0.55, color=bar_colors, edgecolor='gray',
           yerr=bar_err, capsize=3, error_kw={'elinewidth': 1, 'ecolor': 'k'},
           zorder=1)

    rng = np.random.default_rng(seed=0)
    for name in unique_names:
        gdf = element_to_guide_dict.get(name)
        if gdf is None or 'CohesinDependence.perGuide' not in gdf.columns:
            continue
        vals = gdf['CohesinDependence.perGuide'].dropna().values
        if len(vals) == 0:
            continue
        x_jit = positions[name] + rng.uniform(-0.15, 0.15, size=len(vals))
        ax.plot(x_jit, vals, 'o', markersize=dot_size, alpha=dot_alpha,
                color=dot_color, markeredgewidth=0, zorder=2)

    if 'PowerCategory' in gene_enh.columns:
        present_cats = [c for c in gene_enh['PowerCategory'].unique()
                        if c in cohesin_dep_colors]
        handles = [plt.Rectangle((0, 0), 1, 1, color=cohesin_dep_colors[c], ec='gray')
                   for c in present_cats]
        if handles:
            ax.legend(handles, present_cats, frameon=False, fontsize=9, loc='best')

    ax.set_xticks(bar_x)
    if 'DistanceToTSS.Kb' in gene_enh.columns:
        ax.set_xticklabels([f"{d:.1f}" for d in gene_enh['DistanceToTSS.Kb'].values],
                           rotation=0, ha='center')
        ax.set_xlabel('Distance to TSS (Kb)')
    ax.set_ylabel('Cohesin dependence')
    ax.set_title(f'{gene} (n={len(unique_names)} enhancers)')
    ax.axhline(0, linestyle='--', color='gray', lw=1)
    if ylim:
        ax.set_ylim(*ylim)
    ax.spines[['right', 'top']].set_visible(False)

    plt.close()
    return fig, ax



# ------------------ Main ------------------

def main(args):

    # Load data
    print("Loading element list...")
    enhancers= pd.read_csv(args.enhancers, sep=None, engine="python")
    # Add chr, start, end for hg38
    enhancers[['chr', 'start', 'end']] = enhancers['name'].str.split('[:-]', expand=True)
    # Convert start and end to nullable integers
    enhancers['start'] = enhancers['start'].astype('Int64')
    enhancers['end'] = enhancers['end'].astype('Int64')

    print("Loading FF guide levele data...")
    df_guides_noAux = load_df_guides_noAux(args.processed_FF_directory)
    df_guides_plusAux = load_df_guides_plusAux(args.processed_FF_directory)

    # Create element to guide dictionary
    element_to_guide_dict = create_enhancer_to_guide_dict(enhancers, df_guides_noAux, df_guides_plusAux)


    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)


    # Define pilot genes
    pilot_genes = ('CCND1', 'FAM3C', 'MYC', 'SSFA2', 'KITLG')

    # Per-gene TSS midpoint (hg19) for ordering enhancer bars upstream -> TSS -> downstream
    gene_tss_midpoint = {
        'CCND1': (69455280 + 69456448) // 2,
        'FAM3C': (121036122 + 121036876) // 2,
        'KITLG': (88973592 + 88974527) // 2,
        'MYC':   (128745963 + 128748688) // 2,
        'SSFA2': (182755967 + 182757823) // 2,
    }
    # Paper palette from _06 (Fig 4a caption)
    cohesin_dep_colors = {
        'cohesin-dependent': 'darkgreen',
        'cohesin-independent': 'tan',
        'underpowered': 'ivory',
    }

    #Plot Enhancer & Guide Level effects for each gene
    for gene in pilot_genes:
        fig, ax = PlotGuideEffects(enhancers,element_to_guide_dict, gene)
        #Save figure
        path_plot = os.path.join(output_dir, f"CRUDO_{gene}_Enhancer_GuideLevel_Effects.pdf")
        fig.savefig(path_plot, format="pdf", bbox_inches="tight")
        print(f" {gene} plot saved to {path_plot}")

        # v2: element-level CohesinDependence bars with per-gRNA dot overlay
        fig2, ax2 = PlotCohesinDepWithGuides(
            enhancers, element_to_guide_dict, gene,
            cohesin_dep_colors=cohesin_dep_colors,
            gene_tss_midpoint=gene_tss_midpoint,
        )
        path_plot2 = os.path.join(output_dir, f"CRUDO_{gene}_Enhancer_CohesinDep_GuideOverlay.pdf")
        fig2.savefig(path_plot2, format="pdf", bbox_inches="tight")
        print(f" {gene} cohesin-dep overlay saved to {path_plot2}")

if __name__ == "__main__":
    args=parse_args()
    main(args)
