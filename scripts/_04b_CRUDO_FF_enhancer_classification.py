#!/usr/bin/env python3
"""
04b_CRUDO_FF_enhancer_classification.py
Author: Philine Guckelberger
Date: 2025/07/25

Description:
    Script to classify enhancer elements based on statistical significance,
    cohesin dependence, etc. Generates filtered outputs and summary statistics.

Inputs:
    Merged CSV output from 04a_merge_enhancer_gene_datasets.py 
    (equivalent to SupplementaryTable2c).

Outputs:
    - A CSV file consisting of only elements classified as enhancers and their cohesin dependence class
    - A filtered version of this CSV excluding TSS elements
    - Printed summary statistics

Usage:
    python scripts/_04b_CRUDO_FF_enhancer_classification.py \
        --input_file resources/CRUDO_FF_SupplementaryTable2c.csv \
        --output_directory path/to/output/directory/

"""

# ------------------ Imports ------------------

import argparse
import os
import pandas as pd
import numpy as np
import math
import scipy
import statsmodels.stats.multitest as smm
from scipy.stats.contingency import odds_ratio
from scipy.stats import fisher_exact

import warnings
warnings.simplefilter(action="ignore", category=pd.errors.SettingWithCopyWarning)

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Process enhancer data and save results.")
    parser.add_argument("--input_file", help="Path to input CSV file with integrated Hi-C counts, element features, and CRUDO effects from _04a_CRUDO_FF_datasets_merge.py (euivalent to SupplementaryTable2c)")
    parser.add_argument("--output_directory", help="Path to directory to save output files")
    return parser.parse_args()

def load_data(filepath):
    df = pd.read_csv(filepath, sep=None, engine="python")
    return df

def filter_significant_elements(df):
    #Filter for significantly elements (enhancers + TSSs)
    df['delta.noAux'] = abs(df['EnhancerEffect.noAux.Rep1'] - df['EnhancerEffect.noAux.Rep2'])
    median_h3k27 = np.median(df['H3K27ac.RPM.values.noAux'].dropna())

    sig_df = df.loc[
            (df['adj.pval.EnhancerEffect.noAux.Rep1'] < 0.05) & (df['adj.pval.EnhancerEffect.noAux.Rep2'] < 0.05) &
         (df['EnhancerEffect.noAux.Rep1'] >= 0.05) & (df['EnhancerEffect.noAux.Rep2'] >= 0.05) &
        (df['H3K27ac.RPM.values.noAux'] > median_h3k27) &
        (df['delta.noAux'] <= 0.05) &
        (df['n'] >= 10) &
        (
        (df['TargetGene'] != 'KITLG') |
        (
            (df['TargetGene'] == 'KITLG') &
            (df['adj.pval.EnhancerEffect.noAux.Rep3'] < 0.05) &
            (df['adj.pval.EnhancerEffect.noAux.Rep4'] < 0.05) &
            (df['EnhancerEffect.noAux.Rep3'] >= 0.05) &
            (df['EnhancerEffect.noAux.Rep4'] >= 0.05)
        )
         )
            ]
    return sig_df


def add_calculated_columns(sig_df):
    #Add columns for HiC, ABC, H3K27ac changes and perform FDR corrections
    # Hi-C changes
    sig_df['normalized.HiC.5Kb.noAux.pseudo'] = sig_df['SCALE.normalized.observed.5Kb.noAux'] + 1
    sig_df['normalized.HiC.5Kb.Aux.pseudo'] = sig_df['SCALE.normalized.observed.5Kb.Aux.scaled'] + 1
    sig_df['normalized.HiC.5Kb.log2FC.(Aux/noAux)'] = np.log2(sig_df['normalized.HiC.5Kb.Aux.pseudo']) - np.log2(sig_df['normalized.HiC.5Kb.noAux.pseudo'])
    hic_FC = sig_df['normalized.HiC.5Kb.Aux.pseudo'] / sig_df['normalized.HiC.5Kb.noAux.pseudo']
    sig_df['%Change.normalized.HiC.5Kb'] = (hic_FC - 1)

    # ABC changes
    ABC_FC = sig_df['ABC.Score.Aux'] / sig_df['ABC.Score.noAux']
    sig_df['%Change.ABC.Score'] = (ABC_FC - 1)

    # H3K27ac changes
    H3K27ac_FC = sig_df['H3K27ac.RPM.values.Aux'] / sig_df['H3K27ac.RPM.values.noAux']
    sig_df['%Change.H3K27ac.RPM.values'] = (H3K27ac_FC - 1)

    # Generate FDR correction alpha for cohesin dependence to be used in _03_CRUDO_FF_analysis_pipeline.py
    sig_df['adj.pval.CohesinDependence.Rep1'] = smm.fdrcorrection(sig_df["pval.CohesinDependence.Rep1"])[1]
    sig_df['adj.pval.CohesinDependence.Rep2'] = smm.fdrcorrection(sig_df["pval.CohesinDependence.Rep2"])[1]
    sig_df['adj.pval.CohesinDependence'] = smm.fdrcorrection(sig_df["pval.CohesinDependence"])[1]

    rejected, _, _, _ = smm.multipletests(list(sig_df["pval.CohesinDependence"]), method='fdr_bh')
    alpha_adj = 0.05 / sum(rejected)
    print('BH corrected alpha for cohesin dependence - use for power in _03_CRUDO_FF_analysis_pipeline.py: ' + str(alpha_adj))
    print("")

    # Cohesin dependence delta
    sig_df['CohesinDependence.Delta'] = abs(sig_df['CohesinDependence.Rep1'] - sig_df['CohesinDependence.Rep2'])
    return sig_df

def classify_elements(sig_df):
    #Classify elements based on significance, cohesin dependence, enhancer effect, and power
    elements = sig_df['name_hg38'].unique()

    # Significant cohesin-dependent enhancer effect
    for e in elements:
        if (sig_df.loc[sig_df['name_hg38'] == e, 'adj.pval.CohesinDependence'].values[0] <= 0.05) & (sig_df.loc[sig_df['name_hg38'] == e, 'CohesinDependence.Delta'].values[0] <= 0.25):
            sig_df.loc[sig_df['name_hg38'] == e, 'Significant.CohesinDependent.EnhancerEffect(noAux.vs.Aux)'] = True
        else:
            sig_df.loc[sig_df['name_hg38'] == e, 'Significant.CohesinDependent.EnhancerEffect(noAux.vs.Aux)'] = False

    # Increased enhancer effect without auxin
    sig_df['Increased.EnhancerEffect.noAux'] = False
    for e in elements:
        effect_no_aux = sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.noAux'].values[0]
        effect_plus_aux = sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.Aux'].values[0]
        significance = sig_df.loc[sig_df['name_hg38'] == e, 'adj.pval.CohesinDependence'].values[0]
        delta = sig_df.loc[sig_df['name_hg38'] == e, 'CohesinDependence.Delta'].values[0]
        sig_df.loc[sig_df['name_hg38'] == e, 'Increased.EnhancerEffect.noAux'] = (effect_no_aux > effect_plus_aux) & (significance <= 0.05) & (delta <= 0.25)

    # Enhancer effect size range
    sig_df['EnhancerEffectRange'] = ""
    for e in elements:
        enhancer_effect = sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.noAux'].values[0]
        if enhancer_effect < 0.15:
            sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffectRange'] = '<15%'
        elif (0.15 < enhancer_effect < 0.35):
            sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffectRange'] = '15-35%'
        elif (0.35 < enhancer_effect < 0.65):
            sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffectRange'] = '35-65%'
        else:
            sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffectRange'] = '>65%'

    # Element category classification
    for e in elements:
        if sig_df.loc[sig_df['name_hg38'] == e, 'category'].values[0] == 'TSS':
            sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'] = 'TSS'
        elif sig_df.loc[sig_df['name_hg38'] == e, 'CTCF.H3K27acLow'].values[0] == True:
            sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'] = 'CTCF'
        elif (sig_df.loc[sig_df['name_hg38'] == e, 'adj.pval.CohesinDependence'].values[0] <= 0.05) & (sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.noAux'].values[0] > sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.Aux'].values[0]) & (sig_df.loc[sig_df['name_hg38'] == e, 'CohesinDependence.Delta'].values[0] <= 0.25):
            sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'] = 'cohesin-dependent'
        elif (sig_df.loc[sig_df['name_hg38'] == e, 'adj.pval.CohesinDependence'].values[0] <= 0.05) & (sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.noAux'].values[0] <= sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.Aux'].values[0]):
            sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'] = 'cohesin-independent'
        else:
            sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'] = 'other enhancer'

    # Power category
    for e in elements:
        if sig_df.loc[sig_df['name_hg38'] == e, 'category'].values[0] == 'TSS':
            sig_df.loc[sig_df['name_hg38'] == e, 'PowerCategory'] = 'TSS'
        elif (sig_df.loc[sig_df['name_hg38'] == e, 'adj.pval.CohesinDependence'].values[0] <= 0.05) & (sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.noAux'].values[0] > sig_df.loc[sig_df['name_hg38'] == e, 'EnhancerEffect.Aux'].values[0]) & (sig_df.loc[sig_df['name_hg38'] == e, 'CohesinDependence.Delta'].values[0] <= 0.25):
            sig_df.loc[sig_df['name_hg38'] == e, 'PowerCategory'] = 'cohesin-dependent'
        elif (sig_df.loc[sig_df['name_hg38'] == e, 'ElementCategory'].values[0] != 'cohesin-dependent') & ((sig_df.loc[sig_df['name_hg38'] == e, 'PowerToDetect.50%.CohesinDependence'].values[0] < 0.9) | (sig_df.loc[sig_df['name_hg38'] == e, 'CohesinDependence.Delta'].values[0] > 0.25)):
            sig_df.loc[sig_df['name_hg38'] == e, 'PowerCategory'] = 'underpowered'
        else:
            sig_df.loc[sig_df['name_hg38'] == e, 'PowerCategory'] = 'cohesin-independent'

    sig_df = sig_df.loc[sig_df['SCALE.normalized.observed.5Kb.noAux'] >= 5]
    return sig_df

def print_summary_statistics(sig_df_noTSS):
    #Print summary stats on enhancer categories and contact changes
    distal = sig_df_noTSS.loc[sig_df_noTSS['DistanceToTSS.Kb'] > 50]
    proximal = sig_df_noTSS.loc[sig_df_noTSS['DistanceToTSS.Kb'] < 50]
    print('# of enhancers: ' + str(len(sig_df_noTSS)))
    print('# of distal enhancers: ' + str(len(distal)))
    print('# of proximal enhancers: ' + str(len(proximal)))
    print('')

    cohesin_dependent = sig_df_noTSS.loc[sig_df_noTSS['PowerCategory']=='cohesin-dependent']
    cohesin_underpowered = sig_df_noTSS.loc[sig_df_noTSS['PowerCategory']=='underpowered']
    cohesin_independent = sig_df_noTSS.loc[sig_df_noTSS['PowerCategory']=='cohesin-independent']

    cohesin_dependent_distal = cohesin_dependent.loc[cohesin_dependent['DistanceToTSS.Kb']>50]
    cohesin_underpowered_distal = cohesin_underpowered.loc[cohesin_underpowered['DistanceToTSS.Kb']>50]
    cohesin_independent_distal = cohesin_independent.loc[cohesin_independent['DistanceToTSS.Kb']>50]

    cohesin_dependent_prox = cohesin_dependent.loc[cohesin_dependent['DistanceToTSS.Kb']<50]
    cohesin_underpowered_prox = cohesin_underpowered.loc[cohesin_underpowered['DistanceToTSS.Kb']<50]
    cohesin_independent_prox = cohesin_independent.loc[cohesin_independent['DistanceToTSS.Kb']<50]

    print('# of cohesin-dependent enhancers: ' + str(len(cohesin_dependent)))
    print('# of cohesin-dependent distal enhancers: ' + str(len(cohesin_dependent_distal)))
    print('# of cohesin-dependent proximal enhancers: ' + str(len(cohesin_dependent_prox)))

    print('# of cohesin-independent enhancers: ' + str(len(cohesin_independent)))
    print('# of cohesin-independent distal enhancers: ' + str(len(cohesin_independent_distal)))
    print('# of cohesin-independent proximal enhancers: ' + str(len(cohesin_independent_prox)))

    print('# of underpowered enhancers: ' + str(len(cohesin_underpowered)))
    print('# of underpowered distal enhancers: ' + str(len(cohesin_underpowered_distal)))
    print('# of underpowered proximal enhancers: ' + str(len(cohesin_underpowered_prox)))
    print('')

    print('cohesin-dependent avg log2FC contact change: ' + str(np.mean((cohesin_dependent['normalized.HiC.5Kb.log2FC.(Aux/noAux)']))))
    print('cohesin-independent avg log2FC contact change: ' + str(np.mean(cohesin_independent['normalized.HiC.5Kb.log2FC.(Aux/noAux)'])))
    print('underpowered diff avg log2FC contact change: ' + str(np.mean(cohesin_underpowered['normalized.HiC.5Kb.log2FC.(Aux/noAux)'])))

    print('')
    print('cohesin-dependent avg contact percent change: ' + str(np.mean(cohesin_dependent['%Change.normalized.HiC.5Kb'])))
    print('cohesin-independent avg contact percent change: ' + str(np.mean(cohesin_independent['%Change.normalized.HiC.5Kb'])))
    print('underpowered diff avg contact percent change: ' + str(np.mean(cohesin_underpowered['%Change.normalized.HiC.5Kb'])))
    print('')

    large_contact_loss = sig_df_noTSS.loc[sig_df_noTSS['normalized.HiC.5Kb.log2FC.(Aux/noAux)']<-1]
    large_contact_loss_dep = large_contact_loss.loc[(large_contact_loss['PowerCategory']=='cohesin-dependent')]
    large_contact_loss_indep = large_contact_loss.loc[(large_contact_loss['PowerCategory']=='cohesin-independent')]
    large_contact_loss_underpowered = large_contact_loss.loc[(large_contact_loss['PowerCategory']=='underpowered')]
    print('number of E-G pairs with a contact FC more than 2-fold: ' + str(len(large_contact_loss)))
    print('number of E-G pairs with a contact FC more than 2-fold that are cohesin dependent: ' + str(len(large_contact_loss_dep)))
    print('number of E-G pairs with a contact FC more than 2-fold that are cohesin INdependent: ' + str(len(large_contact_loss_indep)))
    print('number of E-G pairs with a contact FC more than 2-fold that are cohesin underpowered: ' + str(len(large_contact_loss_underpowered)))
    print('')

    small_contact_change = sig_df_noTSS.loc[sig_df_noTSS['normalized.HiC.5Kb.log2FC.(Aux/noAux)']>-1]
    small_contact_change_dep = small_contact_change.loc[(small_contact_change['PowerCategory']=='cohesin-dependent')]
    small_contact_change_indep = small_contact_change.loc[(small_contact_change['PowerCategory']=='cohesin-independent')]
    small_contact_change_underpowered = small_contact_change.loc[(small_contact_change['PowerCategory']=='underpowered')]
    print('number of E-G pairs with a contact FC LESS than 2-fold: ' + str(len(small_contact_change)))
    print('number of E-G pairs with a contact FC LESS than 2-fold that are cohesin dependent: ' + str(len(small_contact_change_dep)))
    print('number of E-G pairs with a contact FC LESS than 2-fold that are cohesin INdependent: ' + str(len(small_contact_change_indep)))
    print('number of E-G pairs with a contact FC LESS than 2-fold that are cohesin underpowered: ' + str(len(small_contact_change_underpowered)))
    print('')

    # Compute change in contact
    print('distal avg contact percent change: ' + str(np.mean(distal['%Change.normalized.HiC.5Kb'])))
    print('proximal avg contact percent change: ' + str(np.mean(proximal['%Change.normalized.HiC.5Kb'])))
    print('median enhancer effect: ' + str(np.median(sig_df_noTSS['EnhancerEffect.noAux'])))

    print("")

    #number of enhancers changes in effect between two conditions
    sig_df_noTSS['Reduction.EnhancerEffect'] = sig_df_noTSS['EnhancerEffect.Aux'] / sig_df_noTSS['EnhancerEffect.noAux']
    print('number of enhancers >5% decrease in effect: ' + str(len(sig_df_noTSS.loc[sig_df_noTSS['Reduction.EnhancerEffect']<=0.95])))

    print("")

    # Compute odds ratio for cohesin dependence distal vs proximal
    odds_r= odds_ratio([[len(cohesin_dependent_distal), len(cohesin_dependent_prox)], [len(cohesin_independent_distal),len(cohesin_independent_prox)]])
    print('Odds Ratio cohesin dependence distal vs. proximal: ' + str(odds_r))
    # with Haldane-Anscombe correction
    odds_r_corrected= odds_ratio([[len(cohesin_dependent_distal)+1, len(cohesin_dependent_prox)+1], [len(cohesin_independent_distal)+1,len(cohesin_independent_prox)+1]])
    print('Odds Ratio cohesin dependence distal vs. proximal with Haldane-Anscombe correction: ' + str(odds_r_corrected))
    
    fishers = fisher_exact([[len(cohesin_dependent_distal), len(cohesin_dependent_prox)], [len(cohesin_independent_distal),len(cohesin_independent_prox)]], alternative='two-sided')
    print('Fisher exact cohesin dependence distal vs. proximal: ' + str(fishers))

    # Compute odds ratio for change in 3D contact cohesin depdendent vs independent
    odds_r= odds_ratio([[len(large_contact_loss_dep), len(small_contact_change_dep)], [len(large_contact_loss_indep), len(small_contact_change_indep)]])
    print('Odds Ratio >2-fold change cohesin dependent vs. independent: ' + str(odds_r))
    # with Haldane-Anscombe correction
    odds_r_corrected = odds_ratio([[len(large_contact_loss_dep)+1, len(small_contact_change_dep)+1], [len(large_contact_loss_indep)+1, len(small_contact_change_indep)+1]]) 
    print('Odds Ratio >2-fold change cohesin dependent vs. independent with Haldane-Anscombe correction: ' + str(odds_r_corrected))

    fishers = fisher_exact([[len(large_contact_loss_dep), len(small_contact_change_dep)], [len(large_contact_loss_indep), len(small_contact_change_indep)]], alternative='two-sided')
    print('Fisher exact >2-fold change cohesin dependent vs. independent: ' + str(fishers))



# ------------------ Main ------------------

def main(args):
    # Read in data
    df = load_data(args.input_file)

    # Identify and classify signfiicant elements
    sig_df = filter_significant_elements(df)
    sig_df = add_calculated_columns(sig_df)
    sig_df = classify_elements(sig_df)

    # Generate an enhanceer only (no TSS) version of the df
    sig_df_noTSS = sig_df.loc[sig_df['PowerCategory'] != 'TSS']

    # Print the enhancer stats
    print_summary_statistics(sig_df_noTSS)
    print()

    #make output directory ifit does not exist yet
    output_dir = args.output_directory
    os.makedirs(output_dir, exist_ok=True)

    # Save files
    path_sig_df = os.path.join(output_dir, "CRUDO_FF_significant.csv")
    sig_df.to_csv(path_sig_df, index=False)
    print(f"Saved: {path_sig_df}")

    path_sig_df_noTSS = os.path.join(output_dir, "CRUDO_FF_enhancers.csv")
    sig_df_noTSS.to_csv(path_sig_df_noTSS, index=False)
    print(f"Saved: {path_sig_df_noTSS}")

    
# --- Main execution ---
if __name__ == "__main__":
    args=parse_args()
    main(args)



