#!/usr/bin/env python3

"""
_07_CRUDO_Enhancer_CohesinDep_modeling.py
Author: Philine Guckelberger
Date: 2025/08/06

Description:
Fits models for enhancer Cohesin dependence and plots adjusted R^2 with ANOVA p-values.

Inputs:
    - ACSV file containing CRUDO enhancers and extracted feature data (output from _05b_CRUDO_RAD21_analyzsis.py)

Outputs:
    - Model parameters 
    - Statistics for model comparisons
    - Variance Inflation Factor (VIF) values
    - Adjusted R^2 plot PDF

Usage:
    python scripts/_07_CRUDO_Enhancer_CohesinDep_modeling.py \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --output_directory_plot path/to/output/directory/plots/
"""

# ------------------ Imports ------------------
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.outliers_influence import variance_inflation_factor
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Modeling CRUDO enhancer cohesin dependence ")
    parser.add_argument("--enhancers", required=True, help="CRUDO enhancer csv. Output from _05b_CRUDO_RAD21_analyzsis.py")
    parser.add_argument("--output_directory_plot", required=True, help="Output diretory for plots.")
    return parser.parse_args()


def analyze_and_plot_models(data, response_col, contact_col, distance_col, output_dir):
    # Fit models
    model1 = ols(f"{response_col} ~ {contact_col}", data=data).fit()
    model2 = ols(f"{response_col} ~ {distance_col}", data=data).fit()
    model3 = ols(f"{response_col} ~ {contact_col} + {distance_col}", data=data).fit()
    # Print summaries
    print("Model 1 (Contact Only):")
    print(model1.summary())
    print("\nModel 2 (Distance Only):")
    print(model2.summary())
    print("\nModel 3 (Contact + Distance):")
    print(model3.summary())
    # ANOVA comparisons
    anova_contact_vs_combined = sm.stats.anova_lm(model1, model3, typ=1)
    anova_distance_vs_combined = sm.stats.anova_lm(model2, model3, typ=1)
    # Print ANOVA comparisons
    print("\nANOVA: Model 1 (Contact) vs Model 3 (Combined):")
    print(anova_contact_vs_combined)
    print("\nANOVA: Model 2 (Distance) vs Model 3 (Combined):")
    print(anova_distance_vs_combined)
    # Adjusted R^2 values
    print("\nAdjusted R^2 Values:")
    print(f"Model 1 (Contact Only): {model1.rsquared_adj:.4f}")
    print(f"Model 2 (Distance Only): {model2.rsquared_adj:.4f}")
    print(f"Model 3 (Contact + Distance): {model3.rsquared_adj:.4f}")
    # VIF calculation
    X = data[[contact_col, distance_col]]
    X = sm.add_constant(X)
    vif_data = pd.DataFrame({
        "Variable": X.columns,
        "VIF": [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    })
    print("\nVariance Inflation Factor (VIF):")
    print(vif_data)
    # Prepare data for plotting
    adj_R2_values = {
        "Model 1 (contact)": model1.rsquared_adj, 
        "Model 2 (distance)": model2.rsquared_adj, 
        "Model 3 (combined)": model3.rsquared_adj   
        }
    # Pairwise comparisons vs combined (p-values from ANOVA)
    p_values = {
        "Model 1 vs 3": anova_contact_vs_combined['Pr(>F)'][1],
        "Model 2 vs 3": anova_distance_vs_combined['Pr(>F)'][1]
        }
    anova_df = pd.DataFrame(list(adj_R2_values.items()), columns=['Model', 'adj.R^2'])
    #Plot
    fig, ax = plt.subplots(figsize=(15, 5))
    bars = ax.barh(anova_df['Model'], anova_df['adj.R^2'], color="lightgray", edgecolor='black')
    # Annotate p-values for comparisons vs combined (on the first two bars)
    comparisons = ["Model 1 vs 3", "Model 2 vs 3"]
    for i, comp in enumerate(comparisons):
        ax.text(
            anova_df['adj.R^2'][i] + 0.01, i,  # place to the right of the bar
            f"p={p_values[comp]:.5f}",  # 0.xxxxx format
            va='center', ha='left', fontsize=11, color='black'
            )
    ax.invert_yaxis()  # so the combined model1 is at the top
    ax.set_xlabel('Adjusted R^2', fontsize=12)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #Save figure
    path_plot = os.path.join(output_dir, f"CohesinDep_modeling_adjR2.pdf")
    fig.savefig(path_plot, format="pdf", bbox_inches="tight")
    print(f"Plot saved to {path_plot}")


# ------------------ Main ------------------

def main(args):
    enhancers = pd.read_csv(args.enhancers, sep=None, engine="python")

    # Rename columns for OLS compatibility
    enhancers['contact_change']=enhancers['%Change.normalized.HiC.5Kb']
    enhancers['distance_log10']=np.log10(enhancers['DistanceToTSS.Kb'])
    
    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)

    analyze_and_plot_models(
    data=enhancers,
    response_col='CohesinDependence',
    contact_col='contact_change',
    distance_col='distance_log10',
    output_dir=output_dir)
  
if __name__ == "__main__":
    args = parse_args()
    main(args)