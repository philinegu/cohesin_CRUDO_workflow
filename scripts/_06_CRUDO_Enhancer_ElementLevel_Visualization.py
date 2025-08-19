#!/usr/bin/env python3

"""
06_CRUDO_Enhancer_ElementLevel_Visualization.py
Author: Philine Guckelberger
Date: 2025/07/31

Description:
    This script reads CRUDO enhancer data along with associated feature values
    (e.g., cohesin dependence, Hi-C changes, distances to TSS, RAD21 signal, etc.) and performs
    correlation analyses. It generates a series of scatter plots, grouped plots,
    and bar plots to visualize relationships between features at the enhancer level.

Inputs:
    - A CSV file containing CRUDO enhancers and extracted feature data (output from _05b_CRUDO_RAD21_analysis.py)


Outputs:
    - Scatter plots and bar plots illustrating feature correlations (PDF files of all generated plots, saved to the specified output directory.)
    - Printed Pearson correlations

Usage:
    python scripts/_06_CRUDO_Enhancer_ElementLevel_Visualization.py \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --output_directory_plot path/to/output/directory/plots/

"""

# ------------------ Imports ------------------
import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.lines as mlines
from matplotlib.ticker import FormatStrFormatter
import scipy
from scipy import stats
from scipy.optimize import curve_fit
import statsmodels.api as sm
import statsmodels.stats.multitest as smm
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Visualize CRUDO enhancers")
    parser.add_argument("--enhancers", required=True, help="CRUDO enhancer csv. Output from _05b_CRUDO_RAD21_analysis.py")
    parser.add_argument("--output_directory_plot", required=True, help="Output diretory for plots.")
    return parser.parse_args()



def scatter_with_corr(df, x, y, color,  xlim=(-1,1), ylim=(-1,1),
                      xticks=None, yticks=None, alpha=0.5, s=75, title=None, filename=None, output=None):
    """Generic scatter plot with Pearson correlation printing and saving."""
    fig, ax = plt.subplots(figsize=(5,5))
    df.plot(
        x=x, y=y, kind='scatter',
        lw=1, s=s, legend=False, alpha=alpha, figsize=(5, 5), color=color, ax=ax
    )
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
    # Draw diagonal only if both limits are provided and they are exactly the same
    if xlim != "" and ylim != "" and xlim == ylim:
        plt.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], linestyle='--', color='gray')
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if xticks: ax.set_xticks(xticks)
    if yticks: ax.set_yticks(yticks)
    ax.spines[['right', 'top']].set_visible(False)
    # Compute correlation
    r, p = stats.pearsonr(df[x], df[y])
    if title:
        print(f"{title} R = {r:.5f}, P = {p:.3e}")
    if filename:
        path_plot = os.path.join(output, filename)
        plt.savefig(path_plot, format="pdf", bbox_inches="tight")
    plt.close(fig)
    return (fig, ax) 




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
    # Draw diagonal only if both limits are provided and they are exactly the same
    if xlim != "" and ylim != "" and xlim == ylim:
        plt.plot([xlim[0], xlim[1]], [ylim[0], ylim[1]], linestyle='--', color='gray')
    # Clean plot
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    if xticks: ax.set_xticks(xticks)
    if yticks: ax.set_yticks(yticks)
    ax.spines[['right', 'top']].set_visible(False)
    # Compute correlation
    x_vals = df[x].copy()
    y_vals = df[y].copy()
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



def bar_plot(df, x, y, colors, color_col=None, yerr=None, title="", ylim=None, filename=None, output=None):
    """Bar plot with optional error bars and custom colors."""
    color_values = df[color_col].map(colors) if color_col else 'k'
    fig, ax = plt.subplots(figsize=(8, 4)) 
    df.plot(
        x=x, y=y, kind='bar', fontsize=14,
        color=color_values, legend=False, width=0.5, yerr=yerr, edgecolor='gray', ax=ax
    )
    ax.axhline(y=0, color='k', linestyle='--', lw=1)
    ax.set_title(title)
    if ylim: ax.set_ylim(*ylim)
    ax.spines[['top', 'right']].set_visible(False)
    path_plot = os.path.join(output, filename)
    plt.savefig(path_plot, format="pdf", bbox_inches="tight")
    plt.close(fig)




def ElementLevelRepCorr(enhancers, output_dir, CohesinDepColors, gene):
    #No auxin
    scatter_with_corr(
        df=enhancers,
        x='EnhancerEffect.noAux.Rep1', y='EnhancerEffect.noAux.Rep2', 
        color='#6E7CA0', title=f'{gene} noAux EnhancerEffect Reps', 
        filename=f"CRUDO_{gene}_Enhancers.noAux.RepCorrelation.pdf",
        output=output_dir)
    #Plus auxin
    scatter_with_corr(
        df=enhancers,
        x='EnhancerEffect.Aux.Rep1', y='EnhancerEffect.Aux.Rep2', 
        color='#8F3237', title=f'{gene} aux EnhancerEffect Reps', 
        filename=f"CRUDO_{gene}_Enhancers.Aux.RepCorrelation.pdf",
        output=output_dir)
    #Cohesin dependence
    grouped_scatter(df=enhancers,
        x='CohesinDependence.Rep1', y='CohesinDependence.Rep2', 
        group_col='PowerCategory',colors=CohesinDepColors, title=f'{gene} CohesinDep Reps',
        filename=f"CRUDO_{gene}_Enhancers.CohesinDep.RepCorrelation.pdf", output=output_dir)
    




def PerGeneElementLevelRepCorr(genes,enhancers, output_dir, CohesinDepColors):
    for gene in genes:
        subset=enhancers.loc[enhancers['TargetGene']==gene].copy()
        ElementLevelRepCorr(subset, output_dir, CohesinDepColors, gene)




def CorrLolliplot(corr_dict, filename=None, output=None):
    sorted_correlations = dict(sorted(corr_dict.items(), key=lambda item: abs(item[1][0]),  reverse=False))
    x_sorted = list(sorted_correlations.keys())
    y_sorted =[abs(v.statistic) for v in sorted_correlations.values()]  
    fig, ax = plt.subplots(figsize=(5, 5))
    # Plot the stems
    ax.hlines(y=x_sorted, xmin=0, xmax=y_sorted, color='lightgrey')
    # Plot the markers (the "lollipops")
    ax.plot(y_sorted, x_sorted, 'o', color="k")
    # Set limits and labels
    ax.set_xlim(0, 1)
    ax.set_xlabel("Absolute Pearson's R")
    # Clean up spines
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    path_plot = os.path.join(output, filename)
    plt.savefig(path_plot, format="pdf", bbox_inches="tight")
    plt.close(fig)



# Compute and plot correlations between CohesinDependence and other features
def FeatureCorr(enhancers, CohesinDepColors, output_dir, name=None):
    if name is None:
        name = "Enhancers"
    features=('H3K27ac.RPM.values.noAux', 'SCALE.normalized.observed.5Kb.noAux', 'ABC.Score.noAux',
          'EnhancerEffect.noAux', '%Change.H3K27ac.RPM.values', '%Change.normalized.HiC.5Kb',
           '%Change.ABC.Score', 'DistanceToTSS.Kb', 'RAD21.signal.Aux', '%Change.signal.Rad21' )
    for feature in features:
        if feature =='DistanceToTSS.Kb':
            fig, ax = grouped_scatter(enhancers,
            x="DistanceToTSS.Kb", xlim=(1,1000), xscale="log" ,y='CohesinDependence', ylim=(-0.55,1),
            yerr_col="ci95.CohesinDependence",
            group_col="PowerCategory", colors=CohesinDepColors, title=f'{name} (n={len(enhancers)}):CohesinDependence vs. {feature}')
            ax.axhline(y=50, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1)
            fig.savefig(os.path.join(output_dir, f'CohesinDep_{feature}_scatter_{name}.pdf'))
        else:
            grouped_scatter(enhancers,
            x=feature, y='CohesinDependence', ylim=(-0.55,1),xlim=(""),
            yerr_col="ci95.CohesinDependence",
            group_col='PowerCategory',colors=CohesinDepColors, title=f'{name} (n={len(enhancers)}):CohesinDependence vs. {feature}',
            filename=f'CohesinDep_{feature}_scatter_{name}.pdf', output=output_dir)
    # Build a dictionary of {feature: (r, p)}
    corr_dict = {f: stats.pearsonr(enhancers[f], enhancers['CohesinDependence']) for f in features}
    CorrLolliplot(corr_dict, filename=f"CohesinDep_Feature_Corr_CRUDO_{name}.pdf", output=output_dir)

    

def FeatureColorCorr(enhancers, feat_colors, group_col, output_dir, name=None):
    # CohesinDep vs Change HiC
    grouped_scatter(enhancers,
            x="%Change.normalized.HiC.5Kb", xlim=(-1,0),y='CohesinDependence', ylim=(-0.55,1),
            yerr_col="ci95.CohesinDependence",
            group_col=group_col, colors=feat_colors, title=f'Colored by {name}: CohesinDependence vs. Change in HiC',
            filename=f'Per_{name}_CohesinDep_ChangeHiC_scatter.pdf', output=output_dir)
    # CohesinDep vs Distance
    fig, ax = grouped_scatter(enhancers,
            x="DistanceToTSS.Kb", xlim=(1,1000), xscale="log" ,y='CohesinDependence', ylim=(-0.55,1),
            yerr_col="ci95.CohesinDependence",
            group_col=group_col, colors=feat_colors, title=f'Colored by {name}: CohesinDependence vs. Distance')
    ax.axhline(y=50, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1)
    fig.savefig(os.path.join(output_dir, f'Per_{name}_CohesinDep_Distance_scatter.pdf'))
    # Change HiC vs Distance
    fig, ax = grouped_scatter(enhancers,
            x="%Change.normalized.HiC.5Kb", y='DistanceToTSS.Kb', ylim=(1,1000), yscale="log",xlim=(-1,0), 
            group_col=group_col,colors=feat_colors, title=f'Colored by {name}: Distance vs. change in hic')
    ax.axhline(y=50, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1)
    fig.savefig(os.path.join(output_dir, f'Per_{name}_Distance_ChangeHiC_scatter.pdf'))



def power_law_fit(x, a, b):
    return a * x**b


def fit_power_law_model(x,y):
    params, _ = curve_fit(power_law_fit, x, y)
    y_pred = power_law_fit(x, *params)
    residuals = y - y_pred
    # R-squared calculation
    ss_total = np.sum((y - np.mean(y))**2)
    ss_residual = np.sum(residuals**2)
    r_squared = 1 - (ss_residual / ss_total)
    
    return params, r_squared, residuals


def plot_power_law_fit(x,y, params, x_threshold=0.15):
    a, b = params
    y_threshold = power_law_fit(x_threshold, a, b)
    print(f"ABC score corresponding to {x_threshold*100:.0f}% effect size: {y_threshold:.4f}")
    # Prepare fitted curve
    x_fit = np.linspace(min(x), max(x), 200)
    y_fit = power_law_fit(x_fit, a, b)
    # Plot
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(x, y, color='k', s=10, alpha=0.8)
    ax.plot(x_fit, y_fit, color='r', lw=2)
    # Define axes
    ax.axhline(y=y_threshold, color='k', linestyle='--', lw=1)
    ax.axvline(x=x_threshold, color='k', linestyle='--', lw=1)
    ax.set_yscale('log')
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    # Clean plot
    ax.spines[['right', 'top']].set_visible(False)
    ax.set_xlabel("EnhancerEffect.noAux")
    ax.set_ylabel("ABC.Score.noAux")

    return fig, ax


def plot_residuals(x, residuals):
    # Define plot
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x, residuals, color='k', s=10)
    ax.axhline(y=0, color='r', linestyle='--')
    # set axes
    ax.set_xlabel("EnhancerEffect.noAux")
    ax.set_ylabel("Residuals")
    # Clean plot
    ax.spines[['top', 'right']].set_visible(False)

    return fig, ax


def analyze_enhancer_ABC(output_dir, enhancers_df, x_col='EnhancerEffect.noAux', y_col='ABC.Score.noAux', x_threshold=0.15):
    x = enhancers_df[x_col].values
    y = enhancers_df[y_col].values
    params, r_squared, residuals = fit_power_law_model(x, y)
    # Spearman correlation
    r_s, p_s = stats.spearmanr(x, y)
    print(f"Spearman r = {r_s:.5f}, P = {p_s:.3e}")
    print(f"R-squared = {r_squared:.4f}")
    # Make plots
    fig1, ax1 = plot_power_law_fit(x, y, params, x_threshold)
    fig2, ax2 = plot_residuals(x, residuals)
    # Save plots
    fig1.savefig(f"{output_dir}/ABC_EnhancerEffect_power_law_fit.pdf", format="pdf", bbox_inches="tight")
    fig2.savefig(f"{output_dir}/ABC_EnhancerEffect_residuals_plot.pdf", format="pdf", bbox_inches="tight")
    print(f"Saved plots to {output_dir}")
    #Close plots
    plt.close(fig1)
    plt.close(fig2)


# ------------------ Main ------------------

def main(args):
    # Load enhancer metadata
    enhancers = pd.read_csv(args.enhancers, sep=None, engine="python")
    Distal_filtered_enhancers = enhancers.loc[enhancers['DistanceToTSS.Kb'] >=50].copy()
    
    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)

    # Configure genes and CohesinDep colors
    genes= ('FAM3C', 'CCND1', 'KITLG', 'SSFA2', 'MYC')
    CohesinDepColors={'cohesin-dependent':'darkgreen', 'cohesin-independent': 'tan', 'underpowered':'ivory'}
    
    # Plot ElementLevel replicate correlation across all enhancers
        # Fig.4b,c
    ElementLevelRepCorr(enhancers,output_dir,CohesinDepColors, "All")

    # Plot ElementLevel replicate correlation on a per gene basis (like Fig. 3c,e, Extended Data Fig. 5-8c,e )

    PerGeneElementLevelRepCorr(genes,enhancers, output_dir, CohesinDepColors) 

    # Plot CohesinDependence & Chane in HiC for each enhancer on a per gene basis (like Fig. 3f,g, Extended Data Fig. 5-8f,g, Fig. 5e,f)
    for gene in genes:
        subset = enhancers.loc[enhancers['TargetGene'] == gene].copy()
        subset['exp_obs_ratio'] = subset['SCALE.normalized.expected.5Kb.noAux']/subset['SCALE.normalized.observed.5Kb.noAux']
        subset['%Change.exp_obs']=subset['exp_obs_ratio']-1
    
        bar_plot(subset, x='DistanceToTSS.Kb', y='CohesinDependence',
             colors=CohesinDepColors, color_col='PowerCategory',
             yerr='ci95.CohesinDependence', title='Cohesin-dependence',
             ylim=(-0.5,1), filename=f"CRUDO_EnhancerBar_Dependence_{gene}.pdf", output=output_dir)

        bar_plot(subset, 'DistanceToTSS.Kb', '%Change.normalized.HiC.5Kb', 
             colors=CohesinDepColors, color_col='PowerCategory',
             title='% change Hi-C', ylim=(-1,0.5),
             filename=f"CRUDO_EnhancerBar_HiC_change_{gene}.pdf", output=output_dir)
    
        bar_plot(subset, 'DistanceToTSS.Kb', '%Change.exp_obs', 
             colors=CohesinDepColors, color_col='PowerCategory',
             title='% change expected', ylim=(-1,0.5),
             filename=f"CRUDO_EnhancerBar_exp_change_{gene}.pdf", output=output_dir)        


    # Compute and plot correlations between CohesinDependence and features (like Fig. 5a-c)
    FeatureCorr(enhancers, CohesinDepColors, output_dir, name="AllEnhancers")
    # Print variance explained by hic
    hic_correlation = (stats.pearsonr(enhancers['%Change.normalized.HiC.5Kb'],enhancers['CohesinDependence']))[0]
    Cohesin_Dep_rep_correlation = (stats.pearsonr(enhancers['CohesinDependence.Rep1'],enhancers['CohesinDependence.Rep2']))[0]
    var_explained = (hic_correlation**2)/(Cohesin_Dep_rep_correlation**2)
    print(f'CohesinDep variance explained by change in Hi-C (based on replicate concordance) : {var_explained}')
    print("")

    # Compyte and plot correlations between CohesinDependence and features for distal enhancers only (like Extended Data Fig. 12a-c)
    FeatureCorr(Distal_filtered_enhancers, CohesinDepColors, output_dir, name="DistalEnhancersOnly")



    # Compute and plot correlation between change in H3k27ac and change in 3D contact (like Extended Data Fig. 12d)
    # For all enhancers
    grouped_scatter(enhancers,
            x="%Change.normalized.HiC.5Kb", y='%Change.H3K27ac.RPM.values' ,xlim=(-1,0), ylim=(""),
            group_col='PowerCategory',colors=CohesinDepColors, title=f'All enhancers (n={len(enhancers)}): Change in H3K27ac vs. Change in HiC',
            filename="AllEnhancers_ChangeH3K27_ChangeHiC_scatter.pdf", output=output_dir)
    # For distal enhancers only
    fig, ax = grouped_scatter(Distal_filtered_enhancers,
            x="%Change.normalized.HiC.5Kb", y='%Change.H3K27ac.RPM.values', xlim=(-1,0), ylim=(""),
            group_col='PowerCategory',colors=CohesinDepColors, title=f'Distal enhancers only (n={len(Distal_filtered_enhancers)}): Change in H3K27ac vs. Change in HiC',
            filename="DistalEnhancersOnly_ChangeH3K27_ChangeHiC_scatter.pdf", output=output_dir)


    # Compute and plot correlation between linear genomic distance and change in 3D contact (like Fig. 5i, Fig. 12i)
    # For all enhancers
    fig, ax = grouped_scatter(enhancers,
            x="%Change.normalized.HiC.5Kb", y='DistanceToTSS.Kb', ylim=(1,1000), yscale="log",xlim=(-1,0), 
            group_col='PowerCategory',colors=CohesinDepColors, title=f'All enhancers (n={len(enhancers)}): Distance vs. change in hic')
    ax.axhline(y=50, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1)
    fig.savefig(os.path.join(output_dir, "AllEnhancers_Distance_ChangeHiC_scatter.pdf"))
    # For distal enhancers only
    fig, ax = grouped_scatter(Distal_filtered_enhancers,
            x="%Change.normalized.HiC.5Kb", y='DistanceToTSS.Kb', ylim=(1,1000), yscale="log",xlim=(-1,0), 
            group_col='PowerCategory',colors=CohesinDepColors, title=f'Distal enhancers only (n={len(Distal_filtered_enhancers)}): Distance vs. change in hic')
    ax.axhline(y=50, xmin=0, xmax=1, color='k', linestyle='dashed', lw=1)
    fig.savefig(os.path.join(output_dir, "DistalEnhancersOnly_Distance_ChangeHiC_scatter.pdf"))


#Compute and plot proximity vs EnhancerEffect at baseline (like Extended Data Fig. 11a)
#Baseline Hi-C
    grouped_scatter(enhancers,
            x='SCALE.normalized.observed.5Kb.noAux', y='EnhancerEffect.noAux', xlim=(''), ylim=(''),
            group_col='PowerCategory',colors=CohesinDepColors, title='All enhancers: Baseline HiC vs EnhancerEffecr',
            filename="CRUDO_EnhancerEffect_HiC_scatter.pdf", output=output_dir)
    #Distance
    fig, ax = grouped_scatter(enhancers,
            x='DistanceToTSS.Kb', xlim=(1,1000), xscale="log",y="EnhancerEffect.noAux", ylim=(""), 
            group_col='PowerCategory',colors=CohesinDepColors, title='All enhancers: Distance vs EnhancerEffecr')
    ax.axvline(x=50, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1)
    fig.savefig(os.path.join(output_dir, "CRUDO_EnhancerEffect_Distance_scatter"))


    # Additional RAD21 scatter plots (lke Extended Data Fig. 11g)
    rad21_cols=('RAD21.signal.Aux', '%Change.signal.Rad21' )
    for col in rad21_cols:
        # vs Dsitance
        fig, ax = grouped_scatter(enhancers,
            x='DistanceToTSS.Kb', xlim=(1,1000), xscale="log",y=col, ylim=(""), 
            group_col='PowerCategory',colors=CohesinDepColors, title=f'All enhancers: Distance vs {col}')
        ax.axvline(x=50, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1)
        fig.savefig(os.path.join(output_dir, f"CRUDO_Distance_{col}_scatter.pdf"))
        # vs change HiC
        grouped_scatter(enhancers,
            x="%Change.normalized.HiC.5Kb", y=col, xlim=(-1,0),ylim=(""),
            group_col='PowerCategory',colors=CohesinDepColors, title=f'All enhancers: Change in Hi-C vs {col}',
            filename=f'CRUDO_ChangeHiC_{col}_scatter.pdf', output=output_dir)

    # Computing and Plotting CohesinDep, Change in HiC, and Distance correlations colored by features (like Extended Data Fig. 11b-d, Extended Data Fig. 11h-j, Extended Data Fig. 12f-h)
    #Correlation scatter plots colored by gene
    gene_colors= {'CCND1':'#85B09A', 'KITLG': '#8E9089', 'SSFA2': '#F3E600', 'FAM3C': 'purple', 'MYC': '#6E7CA0'} 
    gene_col="TargetGene"
    FeatureColorCorr(enhancers, gene_colors, gene_col, output_dir, name="Gene")
    #Correlation scatter plots colored by EnhancerEffect range
    effect_range_colors= {'<15%':'#85B09A', '15-35%': '#8E9089', '35-65%': '#F3E600', '>65%': 'purple'} 
    effect_range_col="EnhancerEffectRange"
    FeatureColorCorr(enhancers, effect_range_colors, effect_range_col, output_dir, name="EnhancerEffect")
    #Correlation scatter plots colored by CTCF proximity
    ctcf_colors= {True:'#F3E600', False: '#85B09A'}
    ctcf_col="CTCFwithin5Kb"
    FeatureColorCorr(enhancers, ctcf_colors, ctcf_col, output_dir, name="CTCF")

    # Compute and plot change in contact observed over expected
    enhancers['DistanceBin'] = pd.cut(
            enhancers['DistanceToTSS.Kb'], 
            bins=[0, 50, 100, 500, 5000], 
            labels=["0-50 kb", "50-100 kb", "100-500 kb" ,"500-5000 kb"],
            right=False
            )
    DistanceBinColors = {
            "0-50 kb": "#FBDB93",
            "50-100 kb": "#FFC6C6",
            "100-500 kb": "#8A2D3B",
            "500-5000 kb": "k"
            }   
    fig, ax = grouped_scatter(enhancers,
            x='SCALE.normalized.expected.5Kb.Aux.scaled', xlim=(0.1, 1000),xscale="log", y="SCALE.normalized.observed.5Kb.Aux.scaled", ylim=(0.1,1000), yscale="log", 
            group_col='DistanceBin',colors=DistanceBinColors)
    x_vals = np.linspace(0.1, 1000, 1000)
    ax.plot(x_vals, 2 * x_vals, linestyle='--', color='gray')
    ax.plot(x_vals, 0.5 * x_vals, linestyle='--', color='gray')
    fig.savefig(os.path.join(output_dir, "CRUDO_obs_exp_HiC_scatter.pdf"))

    #Compute # of elements with >2fold higher contact than expected based on distance
    within_counts = {}
    outside_counts = {}
    for status, group in enhancers.groupby('DistanceBin',  observed=False):
        x = group['SCALE.normalized.expected.5Kb.Aux.scaled']
        y = group['SCALE.normalized.observed.5Kb.Aux.scaled']
        fold_change = y / x
        within = (fold_change >= 0.5) & (fold_change <= 2)
        within_count = within.sum()
        outside_count = (~within).sum()
        within_counts[status] = within_count
        outside_counts[status] = outside_count
        print(f"{status}: {within_count} within 2-fold, {outside_count} outside 2-fold")

    #Comput and plot ABC score and enhancer effect correlation (like Extended Data Fig. 14a,b)
    analyze_enhancer_ABC(output_dir, enhancers, x_col='EnhancerEffect.noAux', y_col='ABC.Score.noAux', x_threshold=0.15)


if __name__ == "__main__":
    args = parse_args()
    main(args)