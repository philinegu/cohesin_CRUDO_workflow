#!/usr/bin/env python3
"""
01_get_HiC_contacts.py
Author: Philine Guckelberger
Date: 2025/07/23

Description:
    Extract Hi-C interaction counts for defined genomic regions (i.e. CRUDO element-TSS pairs or ABC enhancer-TSS pairs)
    from .hic files (control and treated conditions). Outputs a combined table of interaction counts for each pair of interest.
Inputs:
    - CSV file with genomic regions of interest. Here a list of all CRUDO non-intra-target genic tested elements.
        Required columns (make sure these are in these have the same Genome ID as the hic files, in this case hg38):
            - chr_hg38: Chromosome (ie., 'chr1')
            - start_hg38: Enhancer start position
            - end_hg38: Enhancer end position
            - TargetGeneTSS_hg38: Transcription start site of the target gene
    - .hic files for untreated and treated conditions. (more efficient if downloaded)
    - Parameters for resolution, normalization, signal type, and interaction window.

Outputs:
    - Per-chromosome interaction files (gzipped TSVs) for untreated and treated samples.
    - A final combined CSV file containing element-TSS pairs with:
        - Extracted Hi-C counts for untreated
        - Extracted Hi-C counts for treated
    Can run multiple times for different normalization or signal types. Need to be merged separately.

Parameters:
    ELEMENT_FILE (str):
        Path to the CSV file containing genomic elements of interest
        Required columns: 'chr_hg38', 'start_hg38', 'end_hg38', 'TargetGeneTSS_hg38'.
        Assuming genomebuild is hg38, otherwise need to change in script configuration

    HIC_UNTREATED (str):
        Path to the untreated condition Hi-C .hic file.  In our case: HCT116 RAD21 dTAG control (ENCFF528XGK.hic)

    HIC_TREATED (str):
        Path to the treated condition Hi-C .hic file. In our case: HCT116 RAD21 dTAG 6h auxin-treated (ENCFF317OIA.hic)

    OUTPUT_DIRECTORY (str):
        Directory where output files will be saved. Per-chromosome files
        are saved under subdirectories 'treated/' and 'untreated/'.

    RESOLUTION (int):
        Resolution in base pairs for Hi-C data extraction (e.g., 5000 for 5Kb).

    SIGNAL_TYPE (str):
        Hi-C signal type to extract.
        Options: 'observed', 'expected', 'oe' (observed over expected).

    NORM_METHOD (str):
        Normalization method to apply.
        Options: 'NONE', 'SCALE', 'Coverage', 'INTER_SCALE', etc.

    MIN_WINDOW (int):
        Minimum genomic distance (bp) between interacting bins to include.

    MAX_WINDOW (int):
        Maximum genomic distance (bp) between interacting bins to include.

    CHROMOSOMES (list):
        List of chromosomes (ints and/or strings) to process (e.g., [1,2,...,22,'X']).
    

Usage:
For parameter explanation:
python _01_get_HiC_contacts.py --help

To explore hi-c parameters:
python scripts/_01_get_HiC_contacts.py --print-metadata /path/to/hic_map.hic

To run extraction:
python scripts/_01_get_HiC_contacts.py \
  --element_file resources/TargetList.csv \
  --hic_untreated path/to/ENCFF528XGK.hic \
  --hic_treated path/to/ENCFF317OIA.hic \
  --output_directory path/to/output/directory/ \
  --resolution 5000 \
  --signal_type observed \
  --norm_method SCALE \
  --min_window 0 \
  --max_window 5000000 \
  --chromosomes 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

"""

# ------------------ Imports ------------------

import os
import sys
import argparse
import numpy as np
import pandas as pd
import hicstraw

# ------------------ HELPER FUNCTIONS ------------------

def print_metadata(hic_file_path):
    # need to import hi-c star seperatly because this is can be run independently of the rest
    import hicstraw
    hic_map = hicstraw.HiCFile(hic_file_path)
    print(f"Genome alignment: {hic_map.getGenomeID()}")
    print(f"Available resolutions: {hic_map.getResolutions()}")
    print("Chromosomes and their lengths:")
    for chrom in hic_map.getChromosomes():
        print(f"  {chrom.name}: {chrom.length}")

def parse_args():
    parser = argparse.ArgumentParser(description="Extract Hi-C interactions for enhancer-TSS pairs or print Hi-C file metadata.")
    parser.add_argument("--element_file", help="CSV file with genomic regions of interest. Columns: chr_hg38, start_hg38, end_hg38, TargetGeneTSS_hg38")
    parser.add_argument("--hic_untreated", help="Path to untreated Hi-C .hic file.")
    parser.add_argument("--hic_treated", help="Path to treated Hi-C .hic file.")
    parser.add_argument("--output_directory", help="Path to output directory for results.")
    parser.add_argument("--resolution", type=int, default=5000, help="Resolution for Hi-C extraction (default: 5000).")
    parser.add_argument("--signal_type", choices=["observed", "expected", "oe"], default="observed", help="Signal type (default: observed).")
    parser.add_argument("--norm_method", default="SCALE", help="Normalization method (default: SCALE).")
    parser.add_argument("--min_window", type=int, default=0, help="Minimum genomic distance (default: 0).")
    parser.add_argument("--max_window", type=int, default=5000000, help="Maximum genomic distance in bp (default: 5000000).")
    parser.add_argument("--chromosomes", nargs="+", default=[str(i) for i in range(1, 23)] + ["X"], help="List of chromosomes (default: 1-22, X).")
    parser.add_argument("--print_metadata", metavar="HIC_FILE", help="Print metadata for a given .hic file and exit.")
    return parser.parse_args()

def save_to_dataframe(result, output):
    binX, binY, counts = [], [], []
    for r in result:
        binX.append(r.binX)
        binY.append(r.binY)
        counts.append(r.counts)
    df = pd.DataFrame({0: binX, 1: binY, 2: counts})
    df.to_csv(output, sep="\t", index=False, header=False)

def load_hic_data(hic_file, chromosome, output_path, args):
    result = hicstraw.straw(args.signal_type, args.norm_method, hic_file, f'chr{chromosome}', f'chr{chromosome}', 'BP', args.resolution)
    save_to_dataframe(result, output_path)
    df = pd.read_csv(output_path, names=['binX', 'binY', 'hic_extracted_counts'], header=None, sep="\t")
    mask = np.logical_and(abs(df['binX'] - df['binY']) <= args.max_window, abs(df['binX'] - df['binY']) >= args.min_window)
    df = df.loc[mask]
    df[['binX', 'binY']] = (df[['binX', 'binY']] / args.resolution).astype(int)
    df['hic_extracted_counts'] = df['hic_extracted_counts'].fillna(0)
    return df

def process_chromosome(chromosome, df_elements, args):
    print(f"Processing chromosome {chromosome}...")
    untreated_path = os.path.join(args.output_directory, 'untreated', f'ENCFF528XGK_chr{chromosome}.{args.resolution}.{args.signal_type}.{args.norm_method}.gz')
    treated_path = os.path.join(args.output_directory, 'treated', f'ENCFF317OIA_chr{chromosome}.{args.resolution}.{args.signal_type}.{args.norm_method}.gz')

    os.makedirs(os.path.dirname(untreated_path), exist_ok=True)
    os.makedirs(os.path.dirname(treated_path), exist_ok=True)

    hic_df = load_hic_data(args.hic_untreated, chromosome, untreated_path, args)
    hic_df_aux = load_hic_data(args.hic_treated, chromosome, treated_path, args)

    #CONFIGURATION
    CHR_COL = 'chr_hg38' # Change if column names in your element file differ
    START_COL = 'start_hg38'
    END_COL = 'end_hg38'
    TSS_COL = 'TargetGeneTSS'  

    temp_df = df_elements.loc[df_elements[CHR_COL] == f'chr{chromosome}'].copy()
    temp_df['bin_TSS'] = np.floor(temp_df[TSS_COL] / args.resolution).astype(int)
    temp_df['enhancer_midpoint'] = (temp_df[START_COL] + temp_df[END_COL]) / 2
    temp_df['bin_enhancer'] = np.floor(temp_df['enhancer_midpoint'] / args.resolution).astype(int)
    temp_df['binX'] = np.amin(temp_df[['bin_TSS', 'bin_enhancer']], axis=1)
    temp_df['binY'] = np.amax(temp_df[['bin_TSS', 'bin_enhancer']], axis=1)

    temp_df = temp_df.merge(hic_df, how='left', on=['binX', 'binY'])
    temp_df = temp_df.merge(hic_df_aux.rename(columns={'hic_extracted_counts': 'hic_extracted_counts_aux'}), how='left', on=['binX', 'binY'])
    
    return temp_df



# ------------------ Main ------------------

def main(args):
    
    # If user only wants metadata
    if args.print_metadata:
        print_metadata(args.print_metadata)
        sys.exit(0)

    # Otherwise, enforce that extraction parameters are set
    if not (args.element_file and args.hic_untreated and args.hic_treated and args.output_directory):
        sys.exit("Error: --element_file, --hic_untreated, --hic_treated, and --output_directory are required when extracting Hi-C contacts.")

    print("Loading element list...")
    df_elements = pd.read_csv(args.element_file, sep=None, engine="python")

    #make output directory if it does not exist yet
    os.makedirs(args.output_directory, exist_ok=True)
    
    print("Running Hi-C extraction...")
    results_list = []
    for chromosome in args.chromosomes:
        results_list.append(process_chromosome(chromosome, df_elements, args))

    df_elements_hic = pd.concat(results_list, ignore_index=True)
    
    # Create dynamic column and output file names based on parameters
    res_kb = args.resolution // 1000
    col_noAux = f"{args.norm_method}.normalized.{args.signal_type}.{res_kb}Kb.noAux"
    col_Aux = f"{args.norm_method}.normalized.{args.signal_type}.{res_kb}Kb.Aux"

    df_elements_hic[col_noAux] = df_elements_hic['hic_extracted_counts'].fillna(0)
    df_elements_hic[col_Aux] = df_elements_hic['hic_extracted_counts_aux'].fillna(0)

    output_file = os.path.join(args.output_directory, f"TargetList.{args.norm_method}.normalized.{args.signal_type}.{res_kb}Kb.csv")
    
    # Save results
    df_elements_hic.to_csv(output_file, header=True, index=False)
    print(f"Saving final file to {output_file}")

if __name__ == "__main__":
    args = parse_args()
    main(args)
