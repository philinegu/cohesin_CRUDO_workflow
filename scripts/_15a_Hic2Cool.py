#!/usr/bin/env python3

"""
15a_Hic2Cool.py
Author: Philine Guckelberger
Date: 2025/08/19

Description:
    Converts a .hic Hi-C file into a .cool file using specified resolution, 
    normalization method, and signal type.

Input:
    A .hic file
    And specify resolution, norm-method, and signal_type

Outpus.
    A .cool file


Usage:
    python scripts/_15a_Hic2Cool.py \
        --hic path/to/.hic \
        --output_directory path/to/output/directory/ \
        --resolution 5000 \
        --signal_type observed \
        --norm_method SCALE

"""

# ------------------ Imports ------------------
import argparse
import os
from pathlib import Path
import pandas as pd
import numpy as np
import hicstraw

# ------------------ HELPER FUNCTIONS ------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Convert .hic to .cool")
    parser.add_argument("--hic", help="Path to .hic file.")
    parser.add_argument("--output_directory", help="Path to output directory for results.")
    parser.add_argument("--resolution", type=int, default=5000, help="Resolution for Hi-C .cool (default: 5000).")
    parser.add_argument("--signal_type", choices=["observed", "expected", "oe"], default="observed", help="Signal type (default: observed).")
    parser.add_argument("--norm_method", default="SCALE", help="Normalization method (default: SCALE).")
    return parser.parse_args()


def ConvertHic2Cool(hic, hic_file, resolution, signal_type, norm_method, cool_file):
    chrom_sizes = pd.Series({chrom.name: chrom.length for chrom in hic.getChromosomes() if chrom.name != "All"})

    # First write the chromosome sizes:
    with open(hic.getGenomeID() + '.size', 'w') as fsize:
        for chrom in hic.getChromosomes():
            if chrom.name != "All":
                fsize.write(f"{chrom.name}\t{chrom.length}\n")
    # Then write the counts in text file:
    with open(cool_file.replace('.cool', ".txt"), 'w') as fo:
        for i in range(len(chrom_sizes)):
            for j in range(i, len(chrom_sizes)):
                chrom1 = chrom_sizes.index[i]
                chrom2 = chrom_sizes.index[j]
                result = hicstraw.straw(signal_type, norm_method, hic_file, chrom1, chrom2, 'BP', resolution)
                for k in range(len(result)):
                    start1 = result[k].binX
                    start2 = result[k].binY
                    value = result[k].counts
                    fo.write(f"{chrom1}\t{start1}\t{start1}\t{chrom2}\t{start2}\t{start2}\t{value:.6f}\n")

    os.system(f"cooler load --count-as-float -f bg2 {hic.getGenomeID()}.size:{resolution} {cool_file.replace('.cool', '.txt')} {cool_file}")

# ------------------ MAIN ------------------

def main(args):
    # Load hic file
    hic_file = args.hic
    hic = hicstraw.HiCFile(args.hic)
    
    # Check if chosen resolution is possible, if not print error
    resolution= args.resolution
    assert resolution in hic.getResolutions(), \
    f"{resolution} is not part of the possible resolutions {','.join(hic.getResolutions())}"

    # Make output directory
    # Set oupput directory
    outdir=args.output_directory
    os.makedirs(outdir, exist_ok=True)

    # Build output cool filename
    hic_basename = Path(args.hic).stem
    cool_file = f"{outdir}/{hic_basename}_{resolution}_{args.norm_method}_{args.signal_type}.cool"

    # Convert to .cool
    ConvertHic2Cool(hic, hic_file, resolution, args.signal_type, args.norm_method, cool_file)
    print(f"Saved .cool to {cool_file}")

if __name__ == "__main__":
    args = parse_args()
    main(args)
