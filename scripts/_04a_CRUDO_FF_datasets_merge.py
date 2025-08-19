#!/usr/bin/env python3
"""
04a_CRUDO_FF_datasets_merge.py
Author: Philine Guckelberger
Date: 2025/07/24

Description:
    This script merges element-TSS data from multiple sources:
        1) Element=TSS Hi-C contact data (observed, expected, OE)
        2) Element feature annotations
        3) CRUDO element effects

Notes:
    This script also rescales the Hi-C counts for the Aux condition based on 
    total read depth to allow direct comparison between untreated (noAux) and 
    treated (Aux) conditions. If using other Hi-C datasets this needs to be adjusted in the script.

Inputs:
    - Merged CSVs obtained from _01_get_HiC_contacts.py containing element-TSS pairs with:
        -  observed, expected, and observed/expected (oe) extracted Hi-C counts at 5Kb resolution for untreated (noAux) and treated (Aux) 
    - Combined CSV file with element features from _02_get_enhancer_features.py
    - A combined CSV from _03_CRUDO_FF_analysis_pipeline.py with CRUDO element effects, cohesin dependence, confidence intervals, and statistical significance, etc.

Outputs:
    - A final merged CSV file with integrated dataset containing element-TSS pairs with Hi-C, feature, and 
    functional data (this corresponds to SupplementaryTable2c in our case)

Usage:

    python scripts/_04a_CRUDO_FF_datasets_merge.py \
        --hic_path path/to/TargetList.SCALE.normalized.merged.5Kb.csv\
        --features_path path/to/output/directory/CRUDO_element_features.csv \
        --elements_effects_path path/to/output/directory/CRUDO_FF_ElementLevel_analysis.csv\
        --output_directory path/to/output/directory/
"""


# ------------------ Imports ------------------

import pandas as pd
import os
import argparse

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Merge hi-c counts, element features, and CRUDO element effects")
    parser.add_argument('--hic_path', required=True, help='Path to CSV file with hi-c counts (either directly from _01_get_HiC_contacts.py or a merged csv of multiple outputs from _01_get_HiC_contacts.py)')
    parser.add_argument('--features_path', required=True, help='Path to CSV file with hi-c counts from _02_get_element_features.py')
    parser.add_argument('--elements_effects_path', required=True, help='Path to CSV file with hi-c counts from _03_CRUDO_FF_analysis_pipeline.py')
    parser.add_argument('--output_directory', required=True, help='Path to directory to save output CSV')
    return parser.parse_args()

def merge_datasets(hic_path, features_path, elements_path):
    # Read data
    df_hic = pd.read_csv(hic_path, sep=None, engine="python")
    df_features = pd.read_csv(features_path, sep=None, engine="python")
    df_elements = pd.read_csv(elements_path, sep=None, engine="python")

    # Scale Hi-C AUX counts based on overall contacts in each conditon
    total_noAux = 1093601688
    total_Aux = 1527291641
    scaling_factor = total_noAux / total_Aux

    hic_aux_cols = [c for c in df_hic.columns if c.endswith(".Aux") and c.startswith("SCALE.normalized")]
    for col in hic_aux_cols:
        new_col = col + ".scaled"
        df_hic[new_col] = df_hic[col] * scaling_factor

    # Define columns to merge on
    merge_cols = ['name', 'name_hg38', 'TargetGene', 'category']

    # Drop overlapping columns
    cols_to_drop_hic = [c for c in df_hic.columns 
                        if (c in df_elements.columns or c in df_features.columns) and c not in merge_cols]
    cols_to_drop_feat = [c for c in df_features.columns 
                         if (c in df_elements.columns or c in df_hic.columns) and c not in merge_cols]

    df_hic_clean = df_hic.drop(columns=cols_to_drop_hic)
    df_features_clean = df_features.drop(columns=cols_to_drop_feat)

    # Merge
    df_merged = pd.merge(df_elements, df_hic_clean, on=merge_cols, how='left')
    df_merged = pd.merge(df_merged, df_features_clean, on=merge_cols, how='left')

    return df_merged


# ------------------ Main ------------------

def main(args):

    df_merged=merge_datasets(args.hic_path, args.features_path, args.elements_effects_path)

    # Save results
    print(df_merged.head())
    output_dir=args.output_directory
    #make output directory if it does not exist yet
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "CRUDO_FF_SupplementaryTable2c.csv")
    df_merged.to_csv(output_path, header=True, index=False)

    print(f"Final merged dataset saved to: {output_path}")

if __name__ == "__main__":
    args=parse_args()
    main(args)