#!/usr/bin/env python3

"""
_15b_CRUDO_PlotBinsHic.py
Author: Philine Guckelberger
Date: 2025/08/19

Description:
    Plot Enhancer-gene pair Hi-C matrices at a fixed window and color scale.

    Notes:
        - Uses cooler "balanced" matrix if available.
        - Works for intra- and inter-chromosomal pairs.
        - Heatmaps share a robust color scale across all pairs.

Input:
    A hic.cool file (fixed resolution, output from _15a_Hic2Cool.py)
    A CSV file with the following columns:
        - chr_enhancer
        - midpoint_enhancer
        - chr_TSS
        - enhancer_TSS
        - name to save plot using the label

Usage:
  python _15b_CRUDO_PlotBinsHic.py \
        --cool path/to/.cool \
        --pairs CRUDO_EnhancerGenePairs.csv \
        --window-bp 20000 \
        --outdir path/to/output/directory/
"""

# ------------------ Imports ------------------
import argparse
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cooler
plt.rcParams['pdf.fonttype'] = 42

# ------------------ HELPER FUNCTIONS ------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Per-pair, fixed-window Hi-C plots from a .cool file")
    parser.add_argument("--cool", required=True, help="path to .cool file (5kb, 10kb, etc)")
    parser.add_argument("--pairs", required=True, help="CSV with columns: chr_enhancer, midpoint_enhancer,chr_TSS, midpoint_TSS, name")
    parser.add_argument("--window-bp", type=int, default=20000, help="Window size in bp (side length of square)")
    parser.add_argument("--outdir", default="eg_plots", help="Directory to save plots")
    parser.add_argument("--max-pairs", type=int, default=None, help="If you want to limit the number of pairs for testing")
    parser.add_argument("--png", action="store_true", help="If you want to also save as PNGs.")
    parser.add_argument("--quantile", type=float, default=0.95, help="If you want to print vmax for a specific quantile")
    return parser.parse_args()

def ShapeMatrix(arr, target_shape):
    trows, tcols = target_shape
    rows, cols = arr.shape
    top = (trows - rows)//2
    bottom = trows - rows - top
    left = (tcols - cols)//2
    right = tcols - cols - left
    if top < 0 or left < 0:
        rstart = (-top) if top < 0 else 0
        rend = rstart + trows
        cstart = (-left) if left < 0 else 0
        cend = cstart + tcols
        return arr[rstart:rend, cstart:cend]
    return np.pad(arr, ((max(0,top), max(0,bottom)), (max(0,left), max(0,right))), mode="constant", constant_values=np.nan)

def ExtractMatrixCooler(clr, chr_enhancer, midpoint_enhancer, chr_TSS, midpoint_TSS, window_bp):
    binsize = clr.binsize
    half_bins = int(round((window_bp / 2) / binsize))
    side = 2 * half_bins + 1
    bin1 = int(midpoint_enhancer // binsize)
    bin2 = int(midpoint_TSS // binsize)
    start1 = max(0, (bin1 - half_bins) * binsize)
    end1   = (bin1 + half_bins + 1) * binsize
    start2 = max(0, (bin2 - half_bins) * binsize)
    end2   = (bin2 + half_bins + 1) * binsize
    #set enhancers as x and TSS as y
    try:
        sub = clr.matrix(balance=True).fetch((chr_TSS,start2,end2), (chr_enhancer,start1,end1))
    except:
        sub = clr.matrix(balance=False).fetch((chr_TSS,start2,end2), (chr_enhancer,start1,end1))

    return ShapeMatrix(np.asarray(sub), (side,side))

def ComputePercentileVmax(matrices, quantile=95):
    all_vals = np.concatenate([m[~np.isnan(m)].ravel() for m in matrices])
    vmax = np.percentile(all_vals, quantile)
    # Optional: round up nicely
    vmax = np.ceil(vmax * 10) / 10
    return vmax

def PlotColorScale(vmin, vmax, cmap, outdir):
    fig, ax = plt.subplots(figsize=(1.2, 4))
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ax, orientation='vertical', label='Interaction count')
    plt.tight_layout()
    path = os.path.join(outdir, "hic_colorbar.pdf")
    fig.savefig(path, format='PDF')
    plt.close(fig)

def PlotMatrix(sub, out_path, kb_half, binsize, center_pos1=None, center_pos2=None, title=None, vmin=None, vmax=None, cmap=None):
    fig, ax = plt.subplots(figsize=(2.6,2.6), dpi=600)
    n_bins = sub.shape[0]
    offsets = np.linspace(-kb_half, kb_half, n_bins)
    extent = (-kb_half, kb_half, -kb_half, kb_half)
    im = ax.imshow(sub, origin="lower", extent=extent, vmin=vmin, vmax=vmax,
                   interpolation="nearest", aspect="equal", cmap=cmap)
    ax.axhline(0, linewidth=0.5, alpha=0.6)
    ax.axvline(0, linewidth=0.5, alpha=0.6)
    if center_pos1 is not None:
        start1 = (center_pos1 - (kb_half*1000))//binsize*binsize
        end1   = (center_pos1 + (kb_half*1000))//binsize*binsize
        xticks = [-kb_half,0,kb_half]
        xlabels = [
            f"{int(xticks[0])} kb ({int(start1/1000)})",
            f"0 kb ({int(center_pos1/1000)})",
            f"{int(xticks[2])} kb ({int(end1/1000)})"
        ]
        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, fontsize=6)
    if center_pos2 is not None:
        start2 = (center_pos2 - (kb_half*1000))//binsize*binsize
        end2   = (center_pos2 + (kb_half*1000))//binsize*binsize
        yticks = [-kb_half,0,kb_half]
        ylabels = [
            f"{int(yticks[0])} kb ({int(start2/1000)})",
            f"0 kb ({int(center_pos2/1000)})",
            f"{int(yticks[2])} kb ({int(end2/1000)})"
        ]
        ax.set_yticks(yticks)
        ax.set_yticklabels(ylabels, fontsize=6, rotation=30)
    # Set axes label
    ax.set_xlabel("Enhancer offset+genomic pos(kb)", fontsize=7)
    ax.set_ylabel("Gene offset+genomic pos (kb)", fontsize=7)
    # Save fig
    if title:
        ax.set_title(title, fontsize=8, pad=4)
    plt.tight_layout(pad=0.5)
    fig.savefig(out_path, bbox_inches="tight")
    return fig, ax
    plt.close(fig)

# ------------------ MAIN ------------------

def main(args):
    # Load hic cooler file
    clr = cooler.Cooler(args.cool)
    binsize = clr.binsize

    #Load E-G pairs
    df=pd.read_csv(args.pairs)

    # Make output directory
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    #For testing the code you can specify max number of pairs
    if args.max_pairs:
        df = df.head(args.max_pairs)
    hic_matrices = [ExtractMatrixCooler(clr, r.chr_enhancer, int(r.midpoint_enhancer), r.chr_TSS, int(r.midpoint_TSS), args.window_bp) for _,r in df.iterrows()]
    
    #Define window offsets (in Kb)
    kb_half = args.window_bp / 2000.0

    if args.quantile:
        # Compute 95th percentile vmax
        quantile = args.quantile * 100
        vmax_quantile = ComputePercentileVmax(hic_matrices, quantile=quantile)
        print(f'The vmax for the {quantile}th quantile is {vmax_quantile}')

    # Define vmin and vmax
    vmin, vmax = 0, 5

    # Generate color map
    cmap = LinearSegmentedColormap.from_list("hic", ["#ffffff","#ffffb2","#fecc5c","#fd8d3c","#e31a1c","#800026"])
    # Plot color scale
    PlotColorScale(vmin, vmax, cmap, outdir)

    # Plot
    for (row, matrix) in zip(df.itertuples(), hic_matrices):
        base_name = f"{row.chr_enhancer}_{row.midpoint_enhancer}__{row.chr_TSS}_{row.midpoint_TSS}"
        #save as svg because lmshow does weird PDF embedding
        out_pdf = outdir / f"{base_name}.svg"
        PlotMatrix(matrix, out_pdf, kb_half, binsize, center_pos1=row.midpoint_enhancer, center_pos2=row.midpoint_TSS, title=row.name, vmin=vmin, vmax=vmax, cmap=cmap)
        if args.png:
            out_png = outdir / f"{base_name}.png"
            PlotMatrix(matrix, out_png, kb_half, binsize, center_pos1=row.midpoint_enhancer, center_pos2=row.midpoint_TSS, title=row.name, vmin=vmin, vmax=vmax, cmap=cmap)

    print(f"Saved {len(hic_matrices)} plots to {outdir}")

if __name__ == "__main__":
    args = parse_args()
    main(args)
