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
  --abc_file resources/Combined_EP_PilotGenes.csv \
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

def intersect_abc_with_elements_closest(df_elements, df_abc):

    # Filter out chr0
    df_elements = df_elements[df_elements['chr_hg38'] != 'chr0'].copy()

    # NEW: stable numeric IDs for elements (so we don't rely on name/coords)
    df_elements = df_elements.reset_index(drop=True).copy()
    df_elements['elem_id'] = df_elements.index.astype(int)

    # Prepare BED files for elements  (include elem_id as 4th col)
    elements_bed = df_elements[['chr_hg38', 'start_hg38', 'end_hg38', 'elem_id']].dropna(
        subset=['chr_hg38', 'start_hg38', 'end_hg38']
    ).copy()
    elements_bed['start_hg38'] = elements_bed['start_hg38'].astype(int)
    elements_bed['end_hg38']   = elements_bed['end_hg38'].astype(int)
    elements_bed_file = 'elements_temp.bed'
    elements_bed.to_csv(elements_bed_file, sep='\t', header=False, index=False)

    # Prepare BED files for ABC (unchanged)
    abc_bed = df_abc[['chr_ABC', 'start_ABC', 'end_ABC']].copy()
    abc_bed['start_ABC'] = abc_bed['start_ABC'].astype(int)
    abc_bed['end_ABC']   = abc_bed['end_ABC'].astype(int)
    abc_bed_file = 'abc_temp.bed'
    abc_bed.to_csv(abc_bed_file, sep='\t', header=False, index=False)

    # Step 1: Try direct intersection first
    elements_bt = pybedtools.BedTool(elements_bed_file)   # B (has elem_id)
    abc_bt = pybedtools.BedTool(abc_bed_file)             # A

    try:
        # NOTE: wo=True returns A cols + B cols + overlap_length
        intersected = abc_bt.intersect(elements_bt, wo=True, F=0.5).to_dataframe(
            names=[
                'chr_ABC', 'start_ABC', 'end_ABC',      # A (3 cols)
                'chr_hg38', 'start_hg38', 'end_hg38', 'elem_id',  # B (4 cols; includes elem_id)
                'overlap_length'                        # overlap
            ]
        ).drop_duplicates()

        # (Optional) Keep your name_* columns for reporting if you want
        intersected['name_ABC'] = (intersected['chr_ABC'].astype(str) + ':' + 
                                   intersected['start_ABC'].astype(str) + '-' + 
                                   intersected['end_ABC'].astype(str))
        intersected['name_hg38'] = (intersected['chr_hg38'].astype(str) + ':' + 
                                    intersected['start_hg38'].astype(str) + '-' + 
                                    intersected['end_hg38'].astype(str))
        
        print(f"Direct intersections found: {len(intersected)}")
        
    except Exception as e:
        print(f"Warning: Direct intersection failed: {e}")
        intersected = pd.DataFrame(columns=[
            'chr_ABC','start_ABC','end_ABC',
            'chr_hg38','start_hg38','end_hg38','elem_id',
            'overlap_length','name_ABC','name_hg38'
        ])

    # Step 2: Find elements that didn't get matched
    # (Optional) Build names on originals (unchanged)
    df_elements['name_hg38'] = (df_elements['chr_hg38'].astype(str) + ':' + 
                                df_elements['start_hg38'].astype(str) + '-' + 
                                df_elements['end_hg38'].astype(str))
    df_abc['name_ABC'] = (df_abc['chr_ABC'].astype(str) + ':' + 
                          df_abc['start_ABC'].astype(str) + '-' + 
                          df_abc['end_ABC'].astype(str))

    # CHANGED: use elem_id from the SAME subset that went to BED
    matched_ids = set(intersected['elem_id'].astype(int).unique()) if not intersected.empty else set()

    # Only consider elements that actually went into the BED (avoid dropna mismatches)
    valid_elem_ids = set(elements_bed['elem_id'].astype(int).unique())

    # Unmatched = valid elements whose elem_id didn't appear in intersection
    need_closest_ids = valid_elem_ids - matched_ids
    unmatched_elements = df_elements[df_elements['elem_id'].isin(need_closest_ids)].copy()

    print(f"Unmatched elements needing closest match: {len(unmatched_elements)}")
    
    # Step 3: Find closest ABC peak for each unmatched element
    closest_matches = []
    
    for _, element in unmatched_elements.iterrows():
        element_chr = element['chr_hg38']
        element_center = (element['start_hg38'] + element['end_hg38']) / 2
        
        # Find ABC peaks on the same chromosome
        abc_same_chr = df_abc[df_abc['chr_ABC'] == element_chr].copy()
        
        if len(abc_same_chr) == 0:
            print(f"Warning: No ABC peaks on {element_chr} for element {element['name_hg38']}")
            continue
            
        # Calculate distances to all ABC peaks on same chromosome
        abc_same_chr['center'] = (abc_same_chr['start_ABC'] + abc_same_chr['end_ABC']) / 2
        abc_same_chr['distance'] = abs(abc_same_chr['center'] - element_center)
        
        # Find closest ABC peak
        closest_abc = abc_same_chr.loc[abc_same_chr['distance'].idxmin()]
        
        # Add to closest matches
        closest_match = {
            'chr_ABC': closest_abc['chr_ABC'],
            'start_ABC': closest_abc['start_ABC'], 
            'end_ABC': closest_abc['end_ABC'],
            'chr_hg38': element['chr_hg38'],
            'start_hg38': element['start_hg38'],
            'end_hg38': element['end_hg38'],
            'overlap_length': 0,  # No actual overlap for closest matches
            'name_ABC': closest_abc['name_ABC'],
            'name_hg38': element['name_hg38'],
            'match_type': 'closest',
            'distance': closest_abc['distance']
        }
        closest_matches.append(closest_match)
    
    # Combine direct intersections with closest matches
    if closest_matches:
        closest_df = pd.DataFrame(closest_matches)
        intersected['match_type'] = 'intersection'
        intersected['distance'] = 0  # Direct overlaps have 0 distance
        
        all_matches = pd.concat([intersected, closest_df], ignore_index=True)
    else:
        intersected['match_type'] = 'intersection' 
        intersected['distance'] = 0
        all_matches = intersected
    
    print(f"Total matches (intersection + closest): {len(all_matches)}")
    print(f"Match types: {all_matches['match_type'].value_counts().to_dict()}")
    
    # Step 4: Merge with original ABC data to get all ABC scores
    abc_merged = pd.merge(df_abc, all_matches, on=['name_ABC', 'chr_ABC', 'start_ABC', 'end_ABC'])

    # Step 5: Merge back with elements - FIX coordinate formatting first
    # Make sure both have consistent coordinate naming (no .0)
    df_elements['name_hg38_clean'] = (
        df_elements['chr_hg38'].astype(str) + ':' +
        pd.to_numeric(df_elements['start_hg38'], errors='coerce').astype('Int64').astype(str) + '-' +
        pd.to_numeric(df_elements['end_hg38'], errors='coerce').astype('Int64').astype(str)
    )

    all_matches['name_hg38_clean'] = (
    all_matches['chr_hg38'].astype(str) + ':' +
    pd.to_numeric(all_matches['start_hg38'], errors='coerce').astype('Int64').astype(str) + '-' +
    pd.to_numeric(all_matches['end_hg38'], errors='coerce').astype('Int64').astype(str)
    )

    # ADD name_hg38_clean to abc_merged before dropping coordinates
    abc_merged['name_hg38_clean'] = (
    abc_merged['chr_hg38'].astype(str) + ':' +
    pd.to_numeric(abc_merged['start_hg38'], errors='coerce').astype('Int64').astype(str) + '-' +
    pd.to_numeric(abc_merged['end_hg38'], errors='coerce').astype('Int64').astype(str)
    )

    abc_merged_no_coords = abc_merged.drop(columns=['chr_hg38', 'start_hg38', 'end_hg38'], errors='ignore')

    # Use the clean names for merging
    merged = pd.merge(df_elements, abc_merged_no_coords, 
                  left_on=['name_hg38_clean', 'TargetGene'], 
                  right_on=['name_hg38_clean', 'TargetGene'], 
                  how='left')
    

    # Check that every element got matched
    unmatched_final = merged[merged['chr_ABC'].isna()]
    if len(unmatched_final) > 0:
        print(f"Warning: {len(unmatched_final)} elements still unmatched after closest search")
        print("Unmatched elements:")
        # Use available columns for debugging
        available_cols = [col for col in ['name_hg38_clean', 'TargetGene', 'chr_hg38'] if col in unmatched_final.columns]
        print(unmatched_final[available_cols].head())
    
    # Add match quality metrics
    merged['match_distance'] = merged['distance']
    merged['is_direct_intersection'] = merged['match_type'] == 'intersection'

    merged['start_hg38'] = pd.to_numeric(merged['start_hg38'], errors='coerce').astype('Int64')
    merged['end_hg38']   = pd.to_numeric(merged['end_hg38'],   errors='coerce').astype('Int64')

    # Rebuild name_hg38 from chr + int start/end (avoids ".0")
    merged['name_hg38'] = (
        merged['chr_hg38'].astype(str) + ':' +
        merged['start_hg38'].astype(str) + '-' +
        merged['end_hg38'].astype(str)
    )

    print(f"Final merged dataset: {len(merged)} rows")
    print(f"Elements with direct intersections: {sum(merged['is_direct_intersection'])}")
    print(f"Elements matched to closest peak: {sum(~merged['is_direct_intersection'])}")
    
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
    df = intersect_abc_with_elements_closest(df_elements, df_abc)

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
