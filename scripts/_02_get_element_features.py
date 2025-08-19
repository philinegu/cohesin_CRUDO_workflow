#!/usr/bin/env python3

"""
02_get_element_features.py
Author: Philine Guckelberger
Date: 2025/07/23

Description:
    Extract enhancer features such as distances, ABC scores, CTCF proximity,
    and H3K27ac enrichment for a list of element-TSS pairs (i.e. CRUDO elements) .

Inputs:
    - CSV file with genomic regions of interest. In our case all CRUDO non-intra-target-genic tested elements.
        Required columns (make sure these are in these have the same Genome build as ChIP data, in this case hg38):
            - chr_hg38: Chromosome (ie., 'chr1')
            - start_hg38: Enhancer start position
            - end_hg38: Enhancer end position
    - Unthresholded ABC score CSV file.
    - CTCF peaks BED file.
    - Untreated and treated H3K27ac BAM files.
        These need to be indexed before use with: samtools index /path/to/bam_file.bam

Outputs:
    - Combined CSV file with enhancer features.

Usage:
python scripts/_02_get_element_features.py \
  --element_file resources/TargetList.csv \
  --abc_file resources/df_ABC_pilot.csv \
  --ctcf_bed resources/ENCFF072BUT_HCT116_CTCF.bed  \
  --h3k27ac_untreated_bam path/to/SRR6164278_H3K27Ac-untreated.srt.nodup.bam \
  --h3k27ac_treated_bam path/to/SRR6164279_H3K27Ac-treated.srt.nodup.bam \
  --output_directory path/to/output/directory/
"""


# ------------------ Imports ------------------

import os
import argparse
import pandas as pd
import numpy as np
import pybedtools
import pysam

# ------------------ HELPER FUNCTIONS ------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Extract enhancer features")
    parser.add_argument('--element_file', required=True, help='CSV file with enhancer elements')
    parser.add_argument('--abc_file', required=True, help='CSV file with ABC scores (treated and untreated combined)')
    parser.add_argument('--ctcf_bed', required=True, help='BED file with CTCF peaks in untreated condition')
    parser.add_argument('--h3k27ac_untreated_bam', required=True, help='BAM file with H3K27ac reads untreated')
    parser.add_argument('--h3k27ac_treated_bam', required=True, help='BAM file with H3K27ac reads treated')
    parser.add_argument('--output_directory', required=True, help='Path to output directory')
    return parser.parse_args()

def calculate_distance(point1, point2):
    return abs(point1 - point2)

def load_element_list(path):
    df = pd.read_csv(path)
    df['enhancer_midpoint_hg38'] = ((df['start_hg38'] + df['end_hg38']) // 2)
    df['DistanceToTSS.Kb'] = ((df['enhancer_midpoint_hg38'] - df['TargetGeneTSS_hg38']).abs()) / 1000
    return df

def load_abc_scores(path):
    df = pd.read_csv(path)
    df.rename(columns={'chr': 'chr_ABC', 'start': 'start_ABC', 'end': 'end_ABC', 'ABC.Score':'ABC.Score.noAux', 'ABC.Score+aux':'ABC.Score.Aux'}, inplace=True)
    return df[['chr_ABC', 'start_ABC', 'end_ABC', 'ABC.Score.noAux', 'ABC.Score.Aux', 'TargetGene']]

def intersect_abc_with_elements(df_elements, df_abc):
    # Filter out chr0
    df_elements = df_elements[df_elements['chr_hg38'] != 'chr0'].copy()

    # Prepare BED files
    elements_bed = df_elements[['chr_hg38', 'start_hg38', 'end_hg38']].dropna().copy()
    elements_bed['start_hg38'] = elements_bed['start_hg38'].astype(int)
    elements_bed['end_hg38'] = elements_bed['end_hg38'].astype(int)
    elements_bed_file = 'elements_temp.bed'
    elements_bed.to_csv(elements_bed_file, sep='\t', header=False, index=False)

    abc_bed = df_abc[['chr_ABC', 'start_ABC', 'end_ABC']].copy()
    abc_bed['start_ABC'] = abc_bed['start_ABC'].astype(int)
    abc_bed['end_ABC'] = abc_bed['end_ABC'].astype(int)
    abc_bed_file = 'abc_temp.bed'
    abc_bed.to_csv(abc_bed_file, sep='\t', header=False, index=False)

    # Intersect using pybedtools
    elements = pybedtools.BedTool(elements_bed_file)
    abc_peaks = pybedtools.BedTool(abc_bed_file)
    intersected = abc_peaks.intersect(elements, wo=True, F=0.5).to_dataframe(names=[
        'chr_ABC', 'start_ABC', 'end_ABC', 'chr_hg38', 'start_hg38', 'end_hg38', 'overlap_length']).drop_duplicates()

    # Create unique keys for merging
    intersected['name_ABC'] = intersected['chr_ABC'].astype(str) + ':' + intersected['start_ABC'].astype(str) + '-' + intersected['end_ABC'].astype(str)
    intersected['name_hg38'] = intersected['chr_hg38'].astype(str) + ':' + intersected['start_hg38'].astype(str) + '-' + intersected['end_hg38'].astype(str)
    df_abc['name_ABC'] = df_abc['chr_ABC'].astype(str) + ':' + df_abc['start_ABC'].astype(str) + '-' + df_abc['end_ABC'].astype(str)

    # Merge ABC with intersected coordinates
    abc_merged = pd.merge(df_abc, intersected, on=['name_ABC', 'chr_ABC', 'start_ABC', 'end_ABC'])
    # Merge back into element list
    merged = pd.merge(df_elements, abc_merged, on=['name_hg38', 'chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene'], how='left')

    # Reorder columns 
    merged = merged[['TargetGene', 'name', 'chr', 'start', 'end','category', 
                     'name_hg38', 'chr_hg38', 'start_hg38', 'end_hg38', 
                     'TargetGeneTSS_hg38', 'DistanceToTSS.Kb', 'ABC.Score.noAux', 'ABC.Score.Aux']]
    return merged

def analyze_ctcf(elements_bedfile, ctcf_bedfile):
    elements = pybedtools.BedTool(elements_bedfile)
    ctcf_peaks = pybedtools.BedTool(ctcf_bedfile)

    elements_list = [(e.chrom, e.start, e.end, e.fields[3]) for e in elements]
    elements_df = pd.DataFrame(elements_list, columns=['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene'])
    elements_df['midpoint'] = (elements_df['start_hg38'] + elements_df['end_hg38']) / 2

    ctcf_list = [(e.chrom, e.start, e.end) for e in ctcf_peaks]
    ctcf_df = pd.DataFrame(ctcf_list, columns=['chr_hg38', 'start_hg38', 'end_hg38'])
    ctcf_df['midpoint'] = (ctcf_df['start_hg38'] + ctcf_df['end_hg38']) / 2

    min_distances = []
    close_to_ctcf = []

    for idx, elem in elements_df.iterrows():
        distances = ctcf_df[ctcf_df['chr_hg38'] == elem['chr_hg38']]['midpoint'].apply(lambda x: calculate_distance(elem['midpoint'], x))
        min_dist = distances.min() if not distances.empty else np.nan
        min_distances.append(min_dist)
        close_to_ctcf.append(min_dist <= 5000 if not np.isnan(min_dist) else False)

    elements_df['DistanceToCTCF'] = min_distances
    elements_df['CTCFwithin5Kb'] = close_to_ctcf

    # Also get overlaps (boolean) using pybedtools intersect
    ctcf_intersect = elements.intersect(ctcf_peaks, wa=True, c=True).to_dataframe(
        names=['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene', 'overlap_count']
            )
    ctcf_intersect['OverlapsCTCF'] = ctcf_intersect['overlap_count'] > 0
   
    print(elements_df.dtypes)
    print(ctcf_intersect.dtypes)

    merged_ctcf = pd.merge(elements_df.drop(columns='midpoint'), ctcf_intersect[['chr_hg38', 'start_hg38', 'end_hg38', 'OverlapsCTCF']], on=['chr_hg38', 'start_hg38', 'end_hg38'])
    return merged_ctcf

def calculate_h3k27ac_rpm(elements_bedfile, bam_file_untreated, bam_file_treated):
    elements = pybedtools.BedTool(elements_bedfile)
    bam_untreated = pysam.AlignmentFile(bam_file_untreated, 'rb')
    bam_treated = pysam.AlignmentFile(bam_file_treated, 'rb')

    total_reads_untreated = bam_untreated.mapped
    total_reads_treated = bam_treated.mapped

    rpm_data = []
    for region in elements:
        chrom, start, end = region.chrom, region.start, region.end
        target_gene = region.fields[3]
        reads_untreated = bam_untreated.fetch(chrom, start, end)
        reads_treated = bam_treated.fetch(chrom, start, end)

        count_untreated = sum(1 for _ in reads_untreated)
        count_treated = sum(1 for _ in reads_treated)

        rpm_untreated = (count_untreated / total_reads_untreated) * 1e6 if total_reads_untreated > 0 else 0
        rpm_treated = (count_treated / total_reads_treated) * 1e6 if total_reads_treated > 0 else 0

        rpm_data.append((chrom, start, end, target_gene, rpm_untreated, rpm_treated))

    rpm_df = pd.DataFrame(rpm_data, columns=['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene','H3K27ac.RPM.values.noAux', 'H3K27ac.RPM.values.Aux'])
    return rpm_df



# ------------------ Main ------------------

def main(args):
    # Load enhancer elements and compute distances
    df_elements = load_element_list(args.element_file)

    # Load ABC scores
    df_abc = load_abc_scores(args.abc_file)

    # Intersect ABC with elements
    df = intersect_abc_with_elements(df_elements, df_abc)

    # Prepare BED file for elements for downstream analysis
    elements_bedfile = 'elements_for_bed_temp.bed'
    df[['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene']].dropna().astype({'start_hg38': int, 'end_hg38': int}).to_csv(elements_bedfile, sep='\t', header=False, index=False)

    # Analyze CTCF
    df_ctcf = analyze_ctcf(elements_bedfile, args.ctcf_bed)

    # Calculate H3K27ac RPM values
    df_h3k27ac = calculate_h3k27ac_rpm(elements_bedfile, args.h3k27ac_untreated_bam, args.h3k27ac_treated_bam)

    # Merge all features into final df
    df_final = df.merge(df_ctcf, on=['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene'], how='left')
    df_final = df_final.merge(df_h3k27ac, on=['chr_hg38', 'start_hg38', 'end_hg38', 'TargetGene'], how='left')
    df_final.drop_duplicates(inplace=True)
 
    # Identify non-enhancers: close to CTCF and low H3K27ac
    median_h3k27ac = np.median(df_h3k27ac['H3K27ac.RPM.values.noAux'].dropna())
    df_final['CTCF.H3K27acLow'] = ((df_final['OverlapsCTCF']) & (df_final['H3K27ac.RPM.values.noAux'] < median_h3k27ac))

    # Save results
    print(df_final.head())

    output_dir=args.output_directory
    #make output directory if it does not exist yet
    os.makedirs(output_dir, exist_ok=True)

    output_path = os.path.join(output_dir, "CRUDO_element_features.csv")
    df_final.to_csv(output_path, header=True, index=False)

    print(f"Saved final element feature table to {output_path }")

if __name__ == '__main__':
    args = parse_args()
    main(args)
