#!/usr/bin/env python3
"""
03_CRUDO_FF_analysis_pipeline.py
Author: Philine Guckelberger
Date: 2025/07/24

Description:
    CRUDO analysis pipeline to process FF gRNA-level outputs across replicates, merge replicates,
    compute enhancer effects, cohesin dependence, and statistical significance 
    at the element level.

Inputs:
    - TargetList.csv: CSV file with genomic regions of interest 
    (in this case FF outputs are in hg19, so we are using the hg19 cooridantes in our TargetList). 
    Modify column names in the script CONFIGURATION as needed.
        Required columns:
            - name (chr:start-end)
            - chr: Chromosome
            - start: Enhancer start position
            - end: Enhancer end position
            - TargetGene: associated gene
    - Per-guide FF outputs for each replicate (individual files per gene and replicate)

Outputs:
    - A combined dataframe with enhancer effects, cohesin dependence, 
      confidence intervals, and statistical significance.

Usage:

- The alpha differential set is based on running _04b_CRUDO_FF_enhancer_classification.py.
- First, run _03_CRUDO_FF_enhancer_classification.py with the default alpha (0.05) for the cohesin dependence power calculation.
- Then run _04a_CRUDO_FF_datasets_merge.py and _04b_CRUDO_FF_enhancer_classification.py to get the Benjamini-Hochberg corrected alpha for testing cohesin dependence.
- After that, re-run _03_CRUDO_FF_analysis_pipeline.py with the obtained alpha.
- Continue running _04a_CRUDO_FF_datasets_merge.py, _04b_CRUDO_FF_enhancer_classification.py and the rest of the pipeline as needed

python scripts/_03_CRUDO_FF_analysis_pipeline.py  \
    --element_file resources/TargetList.csv \
    --processed_FF_directory resources/byExperimentRep/ \
    --output_directory path/to/output/directory/ \
     --alpha_differential 0.0020833333333333333
     
"""

# ------------------ Imports ------------------

import os
import argparse
import pandas as pd
import numpy as np
import scipy
import statsmodels.api as sm
from scipy import stats
import statsmodels.stats.multitest as smm
import statsmodels.formula.api as smf
from statsmodels.stats.power import TTestPower
from statsmodels.stats.power import TTestIndPower
import statsmodels.stats.power as smp
import math

# ------------------ HELPER FUNCTIONS ------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze CRUDO FF data and save element level data.")
    parser.add_argument("--element_file", help="List of elements to aggregate over. Needs: name, chr, start, end, TargetGene")
    parser.add_argument("--processed_FF_directory", help="Directory that leads to the PerExperimentReps, individual gRNA level data per gene, replicate, and condition")
    parser.add_argument("--output_directory", help="Path to directory to save output CSV")
    parser.add_argument('--alpha_differential', type=float, default=0.05,
                        help='Alpha differential threshold (default=0.05). Can be adapted after running _04b_CRUDO_FF_enhancer_classification.py.')
    return parser.parse_args()

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

def create_element_to_guide_dict(df_TargetList, df_guides_noAux, df_guides_plusAux):
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

        if element not in element_to_guide_dict.keys():
            element_to_guide_dict[element] = df_merge.copy()

    print(element_to_guide_dict[element].columns)
    return element_to_guide_dict

def calculate_FF_ElementLevel(pilot_genes, df_TargetList, element_to_guide_dict, alpha_differential_value):
    alpha_threshold = alpha_differential_value
    print(f"Alpha threshold for cohesin dependence testing is {alpha_threshold}")

    FF_ElementLevel = pd.DataFrame()
    for gene in pilot_genes:
        temp_df = pd.DataFrame()
        df_target_elements = df_TargetList.loc[df_TargetList['TargetGene'] == gene]
        target_elements = df_target_elements['name'].unique()
        temp_df['name'] = target_elements
        for element in target_elements:
            temp_df.loc[temp_df['name'] == element, 'name_hg38'] = df_target_elements.loc[df_target_elements['name'] == element, 'name_hg38'].values[0]
            
            temp_df.loc[temp_df['name'] == element, 'n'] = len(element_to_guide_dict[element])
            temp_df.loc[temp_df['name'] == element, 'TargetGene'] = df_target_elements.loc[df_target_elements['name'] == element, 'TargetGene'].values[0]
            temp_df.loc[temp_df['name'] == element, 'category'] = df_target_elements.loc[df_target_elements['name'] == element, 'category'].values[0]
            
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.noAux.Rep1'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep1_noAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.Aux.Rep1'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep1_plusAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.noAux.Rep2'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep2_noAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.Aux.Rep2'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep2_plusAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.noAux.Rep3'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep3_noAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.Aux.Rep3'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep3_plusAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.noAux.Rep4'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep4_noAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.Aux.Rep4'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_rep4_plusAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.noAux'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_noAux']))
            temp_df.loc[temp_df['name'] == element, 'EnhancerEffect.Aux'] = np.mean(1-(element_to_guide_dict[element]['mleAvg_plusAux']))
            
            stdev_noAux_1 = np.std(element_to_guide_dict[element]['mleAvg_rep1_noAux'], ddof=1)
            stdev_noAux_2 = np.std(element_to_guide_dict[element]['mleAvg_rep2_noAux'], ddof=1)
            stdev_noAux = np.std(element_to_guide_dict[element]['mleAvg_noAux'], ddof=1)
            stdev_plusAux_1 = np.std(element_to_guide_dict[element]['mleAvg_rep1_plusAux'], ddof=1)
            stdev_plusAux_2 = np.std(element_to_guide_dict[element]['mleAvg_rep2_plusAux'], ddof=1)
            stdev_plusAux = np.std(element_to_guide_dict[element]['mleAvg_plusAux'], ddof=1)
            
            sqrt_n = math.sqrt(len(element_to_guide_dict[element]))
            
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.noAux.Rep1'] = 1.96 * (stdev_noAux_1 / sqrt_n)
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.Aux.Rep1'] = 1.96 * (stdev_plusAux_1 / sqrt_n)
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.noAux.Rep2'] = 1.96 * (stdev_noAux_2 / sqrt_n)
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.Aux.Rep2'] = 1.96 * (stdev_plusAux_2 / sqrt_n)
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.noAux'] = 1.96 * (stdev_noAux / sqrt_n)
            temp_df.loc[temp_df['name'] == element, 'ci95.EnhancerEffect.Aux'] = 1.96 * (stdev_plusAux / sqrt_n)
            
            # Cohesin dependence
            # rep1
            avg_delta_1 = np.mean(element_to_guide_dict[element]['Delta_aux_rep1'])
            stdev_delta_aux_1 = np.std(element_to_guide_dict[element]['Delta_aux_rep1'], ddof=1)
            ci95_delta_aux_1 = 1.96 * (stdev_delta_aux_1 / sqrt_n)
            avg_guide_1 = np.mean(1 - element_to_guide_dict[element]['mleAvg_rep1_noAux'])
            
            temp_df.loc[temp_df['name'] == element, 'CohesinDependence.Rep1'] = avg_delta_1 / avg_guide_1
            temp_df.loc[temp_df['name'] == element, 'ci95.CohesinDependence.Rep1'] = abs(ci95_delta_aux_1) / avg_guide_1
            
            # rep2
            avg_delta_2 = np.mean(element_to_guide_dict[element]['Delta_aux_rep2'])
            stdev_delta_aux_2 = np.std(element_to_guide_dict[element]['Delta_aux_rep2'], ddof=1)
            ci95_delta_aux_2 = 1.96 * (stdev_delta_aux_2 / sqrt_n)
            avg_guide_2 = np.mean(1 - element_to_guide_dict[element]['mleAvg_rep2_noAux'])
            
            temp_df.loc[temp_df['name'] == element, 'CohesinDependence.Rep2'] = avg_delta_2 / avg_guide_2
            temp_df.loc[temp_df['name'] == element, 'ci95.CohesinDependence.Rep2'] = abs(ci95_delta_aux_2) / avg_guide_2
            
            # combined
            avg_delta = np.mean(element_to_guide_dict[element]['Delta_aux'])
            stdev_delta_aux = np.std(element_to_guide_dict[element]['Delta_aux'], ddof=1)
            ci95_delta_aux = 1.96 * (stdev_delta_aux / sqrt_n)
            avg_guide = 1 - np.mean(element_to_guide_dict[element]['mleAvg_noAux'])
            
            temp_df.loc[temp_df['name'] == element, 'CohesinDependence'] = avg_delta / avg_guide
            temp_df.loc[temp_df['name'] == element, 'ci95.CohesinDependence'] = abs(ci95_delta_aux) / avg_guide
            
            # Significant effect on gene expression
            # Baseline
            sig_element_rep1 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep1_noAux'], element_to_guide_dict[element]['mleAvg_rep1_noAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.noAux.Rep1'] = sig_element_rep1[1]
            sig_element_rep2 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep2_noAux'], element_to_guide_dict[element]['mleAvg_rep2_noAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.noAux.Rep2'] = sig_element_rep2[1]
            
            sig_element_rep1 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep3_noAux'], element_to_guide_dict[element]['mleAvg_rep3_noAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.noAux.Rep3'] = sig_element_rep1[1]
            sig_element_rep2 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep4_noAux'], element_to_guide_dict[element]['mleAvg_rep4_noAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.noAux.Rep4'] = sig_element_rep2[1]
            
            sig_element = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_noAux'], element_to_guide_dict[element]['mleAvg_noAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.noAux'] = sig_element[1]
            
            # +Auxin
            sig_element_aux_rep1 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep1_plusAux'], element_to_guide_dict[element]['mleAvg_rep1_plusAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.Aux.Rep1'] = sig_element_aux_rep1[1]
            sig_element_aux_rep2 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep2_plusAux'], element_to_guide_dict[element]['mleAvg_rep2_plusAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.Aux.Rep2'] = sig_element_aux_rep2[1]
            sig_element_aux_rep1 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep3_plusAux'], element_to_guide_dict[element]['mleAvg_rep3_plusAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.Aux.Rep3'] = sig_element_aux_rep1[1]
            sig_element_aux_rep2 = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_rep4_plusAux'], element_to_guide_dict[element]['mleAvg_rep4_plusAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.Aux.Rep4'] = sig_element_aux_rep2[1]
            sig_element_aux = stats.ttest_ind(element_to_guide_dict[f'nc_{gene}']['mleAvg_plusAux'], element_to_guide_dict[element]['mleAvg_plusAux'], nan_policy='omit', equal_var=True)
            temp_df.loc[temp_df['name'] == element, 'pval.EnhancerEffect.Aux'] = sig_element_aux[1]
            
            # Significant difference +/- cohesin
            # rep1
            t_test_rep1 = stats.ttest_rel(element_to_guide_dict[element]['mleAvg_rep1_noAux'], element_to_guide_dict[element]['mleAvg_rep1_plusAux'], nan_policy='omit')
            temp_df.loc[temp_df['name'] == element, 'pval.CohesinDependence.Rep1'] = t_test_rep1[1]
            # rep2
            t_test_rep2 = stats.ttest_rel(element_to_guide_dict[element]['mleAvg_rep2_noAux'], element_to_guide_dict[element]['mleAvg_rep2_plusAux'], nan_policy='omit')
            temp_df.loc[temp_df['name'] == element, 'pval.CohesinDependence.Rep2'] = t_test_rep2[1]
            # combined
            t_test = stats.ttest_rel(element_to_guide_dict[element]['mleAvg_noAux'], element_to_guide_dict[element]['mleAvg_plusAux'], nan_policy='omit')
            temp_df.loc[temp_df['name'] == element, 'pval.CohesinDependence'] = t_test[1]
            
            # Power to detect 50% cohesin dependence
            sample_size_element = len(element_to_guide_dict[element])
            diff_effect_cohesin_dep_50 = np.mean(1 - element_to_guide_dict[element]['mleAvg_noAux']) * 0.5
            
            # Alpha from reverse BH on only the enhancer elements that we identify later, when running script for the first time use 0.05 as a default and re-run after identifying significant elements and doing the correction
            alpha_differential = alpha_threshold
            
            # 50% cohesin dependence detection power
            power_differential_set_effect = smp.TTestPower().solve_power(
                effect_size=(diff_effect_cohesin_dep_50 / np.std(element_to_guide_dict[element]['mleAvg_noAux'] - element_to_guide_dict[element]['mleAvg_plusAux'])),
                nobs=sample_size_element,
                alpha=alpha_differential,
                power=None,
                alternative='two-sided'
            )
            
            temp_df.loc[temp_df['name'] == element, 'PowerToDetect.50%.CohesinDependence'] = power_differential_set_effect
        FF_ElementLevel = pd.concat([FF_ElementLevel, temp_df], ignore_index=True)
    
    # Enhancer effects FDR correction
    FF_ElementLevel["adj.pval.EnhancerEffect.noAux.Rep1"] = smm.fdrcorrection(FF_ElementLevel["pval.EnhancerEffect.noAux.Rep1"])[1]
    FF_ElementLevel["adj.pval.EnhancerEffect.noAux.Rep2"] = smm.fdrcorrection(FF_ElementLevel["pval.EnhancerEffect.noAux.Rep2"])[1]
    
    FF_ElementLevel.loc[FF_ElementLevel['TargetGene'] == 'KITLG', 'adj.pval.EnhancerEffect.noAux.Rep3'] = smm.fdrcorrection(
        FF_ElementLevel.loc[FF_ElementLevel['TargetGene'] == 'KITLG', 'pval.EnhancerEffect.noAux.Rep3']
    )[1]
    FF_ElementLevel.loc[FF_ElementLevel['TargetGene'] == 'KITLG', 'adj.pval.EnhancerEffect.noAux.Rep4'] = smm.fdrcorrection(
        FF_ElementLevel.loc[FF_ElementLevel['TargetGene'] == 'KITLG', 'pval.EnhancerEffect.noAux.Rep4']
    )[1]
    
    FF_ElementLevel["adj.pval.EnhancerEffect.noAux"] = smm.fdrcorrection(FF_ElementLevel["pval.EnhancerEffect.noAux"])[1]
    # perform FDR correction
    rejected, pvals_corrected, _, _ = smm.multipletests(list(FF_ElementLevel["pval.EnhancerEffect.noAux"]), method='fdr_bh')
    # get adjusted alpha level
    alpha_adj = 0.05 / sum(rejected)

    return FF_ElementLevel



# ------------------ Main ------------------

def main(args):
    # Load data
    print("Loading element list...")
    df_TargetList= pd.read_csv(args.element_file, sep=None, engine="python")
    
    print("Loading FF guide levele data...")
    df_guides_noAux = load_df_guides_noAux(args.processed_FF_directory)
    df_guides_plusAux = load_df_guides_plusAux(args.processed_FF_directory)
    
    # Create element to guide dictionary
    element_to_guide_dict = create_element_to_guide_dict(df_TargetList, df_guides_noAux, df_guides_plusAux)
    
    # Define pilot genes 
    pilot_genes = ('CCND1', 'FAM3C', 'MYC', 'SSFA2', 'KITLG')
    
    # Calculate CRUDO results at element level
    FF_ElementLevel = calculate_FF_ElementLevel(pilot_genes, df_TargetList, element_to_guide_dict, args.alpha_differential)
    
    # Save results
    print(FF_ElementLevel.head())

    output_dir=args.output_directory
    #make output directory if it does not exist yet
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "CRUDO_FF_ElementLevel_analysis.csv")
    FF_ElementLevel.to_csv(output_path, index=False)

if __name__ == "__main__":
    args=parse_args()
    main(args)
