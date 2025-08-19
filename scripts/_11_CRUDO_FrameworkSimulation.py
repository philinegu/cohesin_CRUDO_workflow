#!/usr/bin/env python3
"""
11_CRUDO_FrameworkSimulation.py
Author: Philine Guckelberger
Date: 2025/07/24

Description:
    Simulates the effect of cohesin removal and enhancer modulation on gene expression over time,
    using a multiple enhancer model. This model incorporates the effects of cohesin degradation and 
    dCas9-KRAB mediated enhancer silencing (CRISPRi) on transcription and
    mRNA stability to predict gene expression dynamics.

Notes:
    - Time is modeled in minutes, spanning 0 to 48 hours with 10-minute intervals.
    - Transcription levels and stable mRNA are calculated using exponential decay based on mRNA half-life.
    - The model assumes additive and multiplicative effects of cohesin and enhancer silencing on transcription.

Key Features:
    - Models transcriptional changes due to cohesin degradation and enhancer activity.
    - Incorporates timing and kinetics of auxin-induced cohesin degradation and dCas9-KRAB induction.
    - Computes stable mRNA levels considering mRNA half-life and degradation.
    - Outputs normalized gene expression metrics suitable for comparison to experimental data.
    - Generates plots illustrating gene expression trajectories and bar charts at specified time points.

Parameters:
    - ce (float): Effect of cohesin removal on gene expression in % (as fractions).
    - e (float): Enhancer effect in % (as fractions).
    - ae (float): Effect of Auxin (reduction of contact) on the enhancer effect in % (as fractions).
    - mhl (int, optional): mRNA half-life in minutes (default: 600).
    - ar (int, optional): Auxin rate ain minutes (when do we start removing cohesin) (default: 30).
    - dr (int, optional): dCas9 induction rate in minutes (time until we start silencing enhancers) (default: 1440).
    - atp (int, optional): Auxin induction time point in minutes (when we start adding auxin) (default: 1080).
    - dtp (int, optional): dCas9 induction time point in minutes (when we start adding dCas9) (default: 0).
    - mtp (int, optional): Measurement time point in minutes (default: 1440).


Outputs:
    - PDF plots visualizing simulated gene expression dynamics and bar charts at the specified measurement time point.


Usage:
    # assumes cohesin has a 40% total effect on gene expression
    # assumes the enhancer has a 25% effect on gene expression
    # assumes cohesin has a 90% effect on the enhancer function

    python scripts/_11_CRUDO_FrameworkSimulation.py \
	    --ce 0.4 \
	    --e 0.25 \
	    --ae 0.9 \
	    --output_directory_plot path/to/output/directory/plots/

"""


# ------------------ Imports ------------------
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import math
plt.rcParams['pdf.fonttype'] = 42 

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Simulate the CRUDO experimental framework")
    parser.add_argument('--ce', required=True, type=float, help='Effect of cohesin removal on gene expression in % (as fractions)')
    parser.add_argument('--e', required=True, type=float, help='Enhancer effect in % (as fractions)')
    parser.add_argument('--ae', required=True, type=float, help='Effect of Auxin (reduction of contact) on the enhancer effect in % (as fractions)')
    parser.add_argument('--mhl',type=int,default=600,  help='mRNA half-life in minutes (default: 600)')
    parser.add_argument('--ar',type=int, default=30, help='Auxin rate ain minutes (when do we start removing cohesin) (default: 30)')
    parser.add_argument('--dr',type=int, default=1440, help='Cas9 induction rate in minutes (time until we start silencing enhancers) (default: 1440)')
    parser.add_argument('--atp',type=int, default=1080, help='Auxin induction time point in minutes (when we start adding auxin) (default: 1080)')
    parser.add_argument('--dtp',type=int, default=0, help='dCas9 induction time point in minutes (when we start adding dCas9) (default: 0)')
    parser.add_argument('--mtp',type=int, default=1440, help='Measurement time point in minutes (default: 1440)')
    parser.add_argument('--output_directory_plot', required=True, help='Directory to save output plots')
    return parser.parse_args()


def multiple_enhancer_simulation(ce, e, ae, mhl=600, ar=30, dr=1440, atp=1080, dtp=0, mtp=1440):

    ### Multiple enhancer model
    # Generate the simulation df
    df = pd.DataFrame()

    # COLUMN - Minutes: have the timeframe go from 0-48h (in minutes)
    df['Minutes']= range(0, 2881, 10)

    # COLUMN - Amount of dCas9 based on induction time and rate (capped between 0 and 1, because we can't have more or less than 0 and 100% of cas9)
    krab_ratio= ((df['Minutes'] - dtp) / dr)
    df['Amount of dCas9-KRAB (% of max)'] = krab_ratio.clip(0, 1)

    # COLUMN - Amount of cohesin degradatio nbased on auxin induction time and rate (capped between 0 and 1, because we can't have more or less than 0 and 100% of cohesin degradation)
    auxin_ratio= ((df['Minutes'] - atp) / ar)
    df['Amount of cohesin degradation (% of max)'] = auxin_ratio.clip(0, 1)

    # COLUMN - Transcription level of gene (in %) in the absence of cohesin
    df['Transcription level of gene minus cohesin'] = 1-(df['Amount of cohesin degradation (% of max)']*ce)

    # COLUMN - Transcription level of gene (in %) with crispri
    cas9_level = 1-df['Amount of dCas9-KRAB (% of max)']
    df['Transcription level of gene plus crispri'] = (1-e)  + (e * cas9_level)

    # COLUMN - Transcription level of gene (in %) in the absence of cohesin with crispri
    degron_enhancer_effect = (1 - (df['Amount of cohesin degradation (% of max)'])* ae)
    df['Transcription level of gene minus cohesin plus crispri']=(df['Transcription level of gene minus cohesin']-(e*(1-((df['Amount of cohesin degradation (% of max)'])* ae)))) + (cas9_level*e* degron_enhancer_effect)
    # COLUMN - Stable mRNA level (in units of txn) in no treatment conditions
    mrna_levels= []
    # Iterate through the rows of the DataFrame
    for i in range(len(df)-1):
        x = df['Minutes'][i]
        y = df['Minutes'][i + 1]
        # Calculate the formula for each pair of consecutive rows
        mrna = 1 / (1 - math.exp(-math.log(2) * (y - x) / mhl))
        #Append each row
        mrna_levels.append(mrna)
    # Add a NAN for the last row
    mrna_levels.append(np.nan)
    # Add values to df
    df['Stable mRNA no treatment']= mrna_levels

    # COLUMN - mRNA degradation & Stable mRNA level  (in units of txn) with crispri
    for i in range(len(df)-1):
        x = df['Minutes'][i]
        y = df['Minutes'][i + 1]
        if i == 0:
            list_mrna=[]
            mrna_0 = df['Transcription level of gene plus crispri'][i] / (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_mrna.append(mrna_0)
            list_mrna.extend([None] * (len(df) - 1))
            df['Stable mRNA plus crispri'] = list_mrna

            list_deg=[]
            deg_0 = mrna_0 * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_deg.append(deg_0)
            list_deg.extend([None] * (len(df) - 1))
            df['mRNA degradation plus crispri'] = list_deg  
        else:
            deg_plus = df['Stable mRNA plus crispri'][i-1] * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            df.loc[df['Minutes']== x, 'mRNA degradation plus crispri'] = deg_plus

            mrna_plus = df['Stable mRNA plus crispri'][i-1] + df['Transcription level of gene plus crispri'][i] - df['mRNA degradation plus crispri'][i]
            df.loc[df['Minutes']== x, 'Stable mRNA plus crispri'] = mrna_plus


    # COLUMN - mRNA degradation & Stable mRNA level (in units of txn) in the absence of cohesin
    for i in range(len(df)-1):
        x = df['Minutes'][i]
        y = df['Minutes'][i + 1]
        if i == 0:
            list_mrna=[]
            mrna_0 = df['Transcription level of gene minus cohesin'][i] / (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_mrna.append(mrna_0)
            list_mrna.extend([None] * (len(df) - 1))
            df['Stable mRNA minus cohesin'] = list_mrna

            list_deg=[]
            deg_0 = mrna_0 * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_deg.append(deg_0)
            list_deg.extend([None] * (len(df) - 1))
            df['mRNA degradation minus cohesin'] = list_deg  
        else:
            deg_plus = df['Stable mRNA minus cohesin'][i-1] * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            df.loc[df['Minutes']== x, 'mRNA degradation minus cohesin'] = deg_plus

            mrna_plus = df['Stable mRNA minus cohesin'][i-1] + df['Transcription level of gene minus cohesin'][i] - df['mRNA degradation minus cohesin'][i]
            df.loc[df['Minutes']== x, 'Stable mRNA minus cohesin'] = mrna_plus

    # COLUMN - mRNA degradation & Stable mRNA level (in units of txn) in the absence of cohesin with crispri
    for i in range(len(df)-1):
        x = df['Minutes'][i]
        y = df['Minutes'][i + 1]
        if i == 0:
            list_mrna=[]
            mrna_0 = df['Transcription level of gene minus cohesin plus crispri'][i] / (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_mrna.append(mrna_0)
            list_mrna.extend([None] * (len(df) - 1))
            df['Stable mRNA minus cohesin plus crispri'] = list_mrna

            list_deg=[]
            deg_0 = mrna_0 * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            list_deg.append(deg_0)
            list_deg.extend([None] * (len(df) - 1))
            df['mRNA degradation minus cohesin plus crispri'] = list_deg
        else:
            deg_plus = df['Stable mRNA minus cohesin plus crispri'][i-1] * (1 - math.exp(-math.log(2) * (y - x) / mhl))
            df.loc[df['Minutes']== x, 'mRNA degradation minus cohesin plus crispri'] = deg_plus
            mrna_plus = df['Stable mRNA minus cohesin plus crispri'][i-1] + df['Transcription level of gene minus cohesin plus crispri'][i] - df['mRNA degradation minus cohesin plus crispri'][i]
            df.loc[df['Minutes']== x, 'Stable mRNA minus cohesin plus crispri'] = mrna_plus
        
    # List of columns to normalize to get from stable mRNA to gene expression and add a "normalized" suffix
    columns_to_normalize = [
    'Stable mRNA no treatment',
    'Stable mRNA minus cohesin',
    'Stable mRNA plus crispri',
    'Stable mRNA minus cohesin plus crispri'
    ]

    # Loop through each column and normalize it
    for column in columns_to_normalize:
        first_value = df[column][0]
        normalized_column = column + ' normalized'
        df[normalized_column] = (df[column] / first_value) * 100

    df['gene expression remaining'] = df['Stable mRNA plus crispri normalized']/df['Stable mRNA no treatment normalized']
    df['gene expression remaining minus cohesin'] = df['Stable mRNA minus cohesin plus crispri normalized']/df['Stable mRNA minus cohesin normalized']
    
    return df


def simulate_RNA_lines (df, dtp, atp, mtp):
    fig, ax = plt.subplots(figsize=(8, 4)) 
    df.plot(x='Minutes', y=['Stable mRNA no treatment normalized', 
    'Stable mRNA minus cohesin normalized', 
    'Stable mRNA plus crispri normalized', 
    'Stable mRNA minus cohesin plus crispri normalized'], kind='line',
    fontsize=16, color=['#090580', '#C08261', '#6E7CA0', '#8F3237'], lw=3,
    ax=ax,
    legend=False)


    ax.set_ylim(0, 105)

    ax.axvline(x=dtp, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1, label='Dox induction (0h)')
    ax.axvline(x=atp, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1, label='Auxin induction (18h)')
    ax.axvline(x=mtp, ymin=0, ymax=1, color='k', linestyle='dashed', lw=1, label='Measurement timepoint (24h)')

    ax.spines[['right', 'top']].set_visible(False)

    return fig,ax


def simulate_RNA_mtp_bar (df, mtp):

    fig, ax = plt.subplots(figsize=(6,4))

     # Filter rows where 'Minutes' == mtp
    temp_df = df.loc[df["Minutes"] == mtp, ['Stable mRNA plus crispri normalized', 'Stable mRNA minus cohesin plus crispri normalized']]
    
    # Plot bar chart of filtered data
    temp_df.plot(kind='bar', color=['#6E7CA0', '#8F3237'], ax=ax)

    nc_noAux = df.loc[df["Minutes"] == mtp, ['Stable mRNA no treatment normalized']].values[0]
    ax.axhline(y=nc_noAux, xmin=0, xmax=1, color='#090580', linestyle='dashed')

    nc_Aux = df.loc[df["Minutes"] == mtp, ['Stable mRNA minus cohesin normalized']].values[0]
    ax.axhline(y=nc_Aux, xmin=0, xmax=1, color='#C08261', linestyle='dashed')

    ax.spines[['right', 'top']].set_visible(False)
    ax.set_ylim(0, 105)

    print(temp_df)
    return fig,ax


def simulate_normalized_RNA_mtp_bar (df, mtp):

    fig, ax = plt.subplots(figsize=(6,4))

     # Filter rows where 'Minutes' == mtp
    temp_df = df.loc[df["Minutes"] == mtp, ['gene expression remaining', 'gene expression remaining minus cohesin']]

    # Plot bar chart of filtered data
    temp_df.plot(kind='bar', color=['#6E7CA0', '#8F3237'], ax=ax)


    ax.axhline(y=1, xmin=0, xmax=1, color='k', linestyle='dashed')

    ax.spines[['right', 'top']].set_visible(False)

    print(temp_df)
    return fig,ax

# ------------------ Main ------------------

def main(args):

    # Set oupput directory
    output_dir=args.output_directory_plot
    os.makedirs(output_dir, exist_ok=True)
    safe_args = "_".join(f"{k}{v}" for k, v in vars(args).items() if k != "output_directory_plot")

    simulation_df = multiple_enhancer_simulation(
        ce=args.ce,
        e=args.e,
        ae=args.ae,
        mhl=args.mhl,
        ar=args.ar,
        dr=args.dr,
        atp=args.atp,
        dtp=args.dtp,
        mtp=args.mtp
         )

    fig, ax = simulate_RNA_lines (df=simulation_df, dtp=args.dtp, atp=args.atp, mtp=args.mtp)
    filename = f"CRUDO_Simulation_Lines_{safe_args}.pdf"
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, filename))

    fig, ax = simulate_RNA_mtp_bar(df=simulation_df, mtp=args.mtp)
    filename = f"CRUDO_Simulation_bar_{safe_args}.pdf"
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, filename))

    fig, ax = simulate_normalized_RNA_mtp_bar(df=simulation_df, mtp=args.mtp)
    filename = f"CRUDO_Simulation_norm_bar_{safe_args}.pdf"
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, filename))
    

if __name__ == "__main__":
    args=parse_args()
    main(args)