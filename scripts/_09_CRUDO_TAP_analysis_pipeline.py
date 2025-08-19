#!/usr/bin/env python3
"""
09_CRUDO_TAP_analysis_pipeline.py
Author: Philine Guckelberger
Date: 2025/07/24

Description:
    CRUDO analysis pipeline to process DC-TAP-seq gRNA-level outputs across replicates, merge replicates,
    compute enhancer effects, cohesin dependence, and statistical significance 
    at the element level.

Inputs:
    - Guide-level processed DC-TAP-seq output file:
        A merged data file containing gRNA effects across multiple replicates and conditions 
        (untreated and auxin-treated), including:
            - gRNA metadata (e.g. sequence, coordinates, target gene)
            - gRNA effect values (per replicate and condition)
            - Element assignment (the enhancer/element each gRNA maps to)
Outputs:
    - A combined dataframe with enhancer effects, cohesin dependence, 
      confidence intervals, and statistical significance.
Usage:
    python scripts/_09_CRUDO_TAP_analysis_pipeline.py  \
        --processed_guides_file resources/CRUDO_tapseq.results.guide_level.txt \
        --output_directory ppath/to/output/directory
"""

# ------------------ Imports ------------------

import os
import argparse
import pandas as pd
import numpy as np
import scipy
import statsmodels.api as sm
from scipy import stats, optimize
import statsmodels.stats.multitest as smm
import statsmodels.formula.api as smf
from statsmodels.stats.power import TTestPower
from statsmodels.stats.power import TTestIndPower
import statsmodels.stats.power as smp
import math

# ------------------ HELPER FUNCTIONS ------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze CRUDO DC-TAP-seq data and save element level data.")
    parser.add_argument("--processed_guides_file", help="DC-TA-seq guide level results CSV")
    parser.add_argument("--output_directory", help="Directory to save output CSV")
    return parser.parse_args()


def calculate_TAP_ElementLevel(df_GuideEffects):
    # Defining elements & genes
    elements= df_GuideEffects['element'].unique()
    genes= df_GuideEffects['gene'].unique()
    TAP_ElementLevel= pd.DataFrame()

    for gene in genes:
        temp_df= pd.DataFrame()
        df_target_elements=df_GuideEffects.loc[df_GuideEffects['gene']==gene]
        target_elements = df_target_elements['element'].unique()
        temp_df['name']= target_elements 
        for element in target_elements:
            target_temp = df_target_elements.loc[df_target_elements['element']==element]
            temp_df[['chr', 'start', 'end']]= temp_df['name'].str.split('[:-]', expand=True)
            temp_df[['TargetGene']]= gene
            temp_df.loc[temp_df['name']== element, 'n']= len(target_temp)

            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.noAux.Rep1'] = np.mean(1-(target_temp['minus_v1_normalized']))
            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.Aux.Rep1'] = np.mean(1-(target_temp['plus_v1_normalized']))
            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.noAux.Rep2'] = np.mean(1-(target_temp['minus_v2_normalized']))
            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.Aux.Rep2'] = np.mean(1-(target_temp['plus_v2_normalized']))
            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.noAux'] = np.mean(1-(target_temp['minus_auxin_effect']))
            temp_df.loc[temp_df['name']== element, 'EnhancerEffect.Aux'] = np.mean(1-(target_temp['plus_auxin_effect']))
        
            stdev_noAux_1 = np.std(target_temp['minus_v1_normalized'], ddof=1)
            stdev_plusAux_1 = np.std(target_temp['plus_v1_normalized'], ddof=1)
            stdev_noAux_2 = np.std(target_temp['minus_v2_normalized'], ddof=1)
            stdev_plusAux_2 = np.std(target_temp['plus_v2_normalized'], ddof=1)
            stdev_noAux =  np.std(target_temp['minus_auxin_effect'], ddof=1)
            stdev_plusAux = np.std(target_temp['plus_auxin_effect'], ddof=1)

            sqrt_n = math.sqrt(len(target_temp))

            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.noAux.Rep1'] = 1.96*(stdev_noAux_1/sqrt_n)
            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.Aux.Rep1'] = 1.96*(stdev_plusAux_1/sqrt_n)
            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.noAux.Rep2'] = 1.96*(stdev_noAux_2/sqrt_n)
            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.Aux.Rep2'] = 1.96*(stdev_plusAux_2/sqrt_n)
            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.noAux'] = 1.96*(stdev_noAux/sqrt_n)
            temp_df.loc[temp_df['name']== element, 'ci95.EnhancerEffect.Aux'] = 1.96*(stdev_plusAux/sqrt_n)

            #Cohesin dependence
            # rep1
            avg_delta_1=np.mean(target_temp['Delta_aux_rep1'])
            stdev_delta_aux_1 = np.std(target_temp['Delta_aux_rep1'],ddof=1)
            ci95_delta_aux_1 = 1.96*(stdev_delta_aux_1/sqrt_n)

            avg_guide_1 = 1-(np.mean(target_temp['minus_v1_normalized']))

            temp_df.loc[temp_df['name']== element, 'CohesinDependence.Rep1'] = avg_delta_1/avg_guide_1
            temp_df.loc[temp_df['name']== element, 'ci95.CohesinDependence.Rep1'] = abs(ci95_delta_aux_1)/avg_guide_1

            # rep2
            avg_delta_2=np.mean(target_temp['Delta_aux_rep2'])
            stdev_delta_aux_2 = np.std(target_temp['Delta_aux_rep2'],ddof=1)
            ci95_delta_aux_2 = 1.96*(stdev_delta_aux_2/sqrt_n)
    
            avg_guide_2 = 1-(np.mean(target_temp['minus_v2_normalized']))

            temp_df.loc[temp_df['name']== element, 'CohesinDependence.Rep2'] = avg_delta_2/avg_guide_2
            temp_df.loc[temp_df['name']== element, 'ci95.CohesinDependence.Rep2'] = abs(ci95_delta_aux_2)/avg_guide_2

            # combined
            avg_delta=np.mean(target_temp['Delta_aux'])
            stdev_delta_aux = np.std(target_temp['Delta_aux'],ddof=1)
            ci95_delta_aux = 1.96*(stdev_delta_aux/sqrt_n)
    
            avg_guide = 1-(np.mean(target_temp['minus_auxin_effect']))

            temp_df.loc[temp_df['name']== element, 'CohesinDependence'] = avg_delta/avg_guide
            temp_df.loc[temp_df['name']== element, 'ci95.CohesinDependence'] = abs(ci95_delta_aux)/avg_guide


            #Significant effect on gene expression
            temp_nc = df_target_elements.loc[df_target_elements['locus'] =='negative_control']
            #Baseline
            sig_element_rep1 = stats.ttest_ind(temp_nc['minus_v1_normalized'], target_temp['minus_v1_normalized'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.noAux.Rep1'] = sig_element_rep1[1]
            sig_element_rep2 = stats.ttest_ind(temp_nc['minus_v2_normalized'], target_temp['minus_v2_normalized'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.noAux.Rep2'] = sig_element_rep2[1]
            sig_element = stats.ttest_ind(temp_nc['minus_auxin_effect'], target_temp['minus_auxin_effect'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.noAux'] = sig_element[1]
            #+Auxin
            sig_element_aux_rep1 = stats.ttest_ind(temp_nc['plus_v1_normalized'], target_temp['plus_v1_normalized'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.Aux.Rep1'] = sig_element_aux_rep1[1]
            sig_element_aux_rep2 = stats.ttest_ind(temp_nc['plus_v2_normalized'], target_temp['plus_v2_normalized'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.Aux.Rep2'] = sig_element_aux_rep2[1]
            sig_element_aux = stats.ttest_ind(temp_nc['plus_auxin_effect'], target_temp['plus_auxin_effect'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name']== element, 'pval.EnhancerEffect.Aux'] = sig_element_aux[1]


            #Significant difference +/-cohesin
            #rep1
            t_test = stats.ttest_rel(target_temp['minus_v1_normalized'], target_temp['plus_v1_normalized'], nan_policy='omit')
            temp_df.loc[temp_df['name']== element, 'pval.CohesinDependence.Rep1'] = t_test[1]
            #rep2
            t_test = stats.ttest_rel(target_temp['minus_v2_normalized'], target_temp['plus_v2_normalized'], nan_policy='omit')
            temp_df.loc[temp_df['name']== element, 'pval.CohesinDependence.Rep2'] = t_test[1]
            #combined
            t_test = stats.ttest_rel(target_temp['minus_auxin_effect'], target_temp['plus_auxin_effect'], nan_policy='omit')
            temp_df.loc[temp_df['name']== element, 'pval.CohesinDependence'] = t_test[1]

        if len(TAP_ElementLevel)<0:
            TAP_ElementLevel=temp_df.copy()
        else:
            TAP_ElementLevel = pd.concat([TAP_ElementLevel, temp_df])

    # Enhancer effects FDR correction
    TAP_ElementLevel["adj.pval.EnhancerEffect.noAux.Rep1"]=smm.fdrcorrection(TAP_ElementLevel["pval.EnhancerEffect.noAux.Rep1"])[1]
    TAP_ElementLevel["adj.pval.EnhancerEffect.noAux.Rep2"]=smm.fdrcorrection(TAP_ElementLevel["pval.EnhancerEffect.noAux.Rep2"])[1]
    TAP_ElementLevel["adj.pval.EnhancerEffect.noAux"] = smm.fdrcorrection(TAP_ElementLevel["pval.EnhancerEffect.noAux"])[1]
    # perform FDR correction
    rejected, pvals_corrected, _, _ = smm.multipletests(list(TAP_ElementLevel["pval.EnhancerEffect.noAux"]), method='fdr_bh')
    # get adjusted alpha level
    alpha_adj = 0.05 / sum(rejected)

    return TAP_ElementLevel


# ------------------ Main ------------------

def main(args):
    # Load data
    print("Loading DC-TAP guide level data...")
    df_GuideEffects= pd.read_csv(args.processed_guides_file, sep=None, engine="python")
    
    # Add delta aux columns
    df_GuideEffects['Delta_aux']= (1-(df_GuideEffects['minus_auxin_effect'])) - (1-(df_GuideEffects['plus_auxin_effect']))
    df_GuideEffects['Delta_aux_rep1']= (1-(df_GuideEffects['minus_v1_normalized'])) - (1-(df_GuideEffects['plus_v1_normalized']))
    df_GuideEffects['Delta_aux_rep2']= (1-(df_GuideEffects['minus_v2_normalized'])) - (1-(df_GuideEffects['plus_v2_normalized']))

    
    # Calculate CRUDO results at element level
    TAP_ElementLevel = calculate_TAP_ElementLevel(df_GuideEffects)
    
    # Print results
    print(TAP_ElementLevel.head())

    #make output directory if it does not exist yet
    output_dir=args.output_directory
    os.makedirs(output_dir, exist_ok=True)

    #save results
    output_path = os.path.join(output_dir, "CRUDO_TAP_ElementLevel_analysis.csv")
    TAP_ElementLevel.to_csv(output_path, index=False)

if __name__ == "__main__":
    args=parse_args()
    main(args)
