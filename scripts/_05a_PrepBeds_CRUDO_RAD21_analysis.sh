#!/bin/bash
"""
05a_PrepBeds_CRUDO_RAD21_analysis.sh
Author: Philine Guckelberger
Date: 2025/07/30

Description:
  Prepares BED files and computes average RAD21 ChIP-seq signal over CRUDO enhancers and 
  CTCF sites within target gene loci for subsequent python anaylyses (05b__CRUDO_RAD21_analysis.py). Steps include:
    - Sorting and reducing CTCF BED to essential columns
    - Liftover of hg38 CTCF BED to hg19
    - Extending TSS coordinates by ±1 Mb
    - Intersecting CTCF sites with extended TSS regions to get target-gene-associated CTCF sites
    - Converting CTCF and enhancer BEDs to 4-column format for bigWigAverageOverBed
    - Computing average BigWig signals over enhancers and CTCF sites for treated and untreated conditions

Inputs:
    - CTCF BED file in hg38 (e.g., ENCFF072BUT_HCT116_CTCF.bed)
    - CRUDO enhancer CSV (e.g., 04b_CRUDO_FF_enhancers.csv)
    - TSS BED file in hg19 (e.g., Pilot_TSS_hg19.txt)
    - RAD21 ChIP-seq bigWig files (treated and untreated)
    - LiftOver chain file (hg38ToHg19.over.chain)
    - bedtools and bigWigAverageOverBed installed and available in PATH



Outputs:
    - Sorted and lifted CTCF BED files (hg19)
    - Extended TSS BED files (±1Mb)
    - Intersected CTCF sites within target gene loci BED
    - 4-column BED files for CTCF and enhancers
    - Average BigWig signal files over enhancers and CTCF sites for each condition
 

Usage:
    scripts/_05a_PrepBeds_CRUDO_RAD21_analysis.sh
"""

# Sort and keep first 3 columns of CTCF BED
awk '{print $1, $2, $3}' resources/ENCFF072BUT_HCT116_CTCF.bed | sort -k1,1 -k2,2n > path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted.bed

#liftover hg38 BED file to hg19
path/to/liftOver \
  path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted.bed \
  path/to/hg38ToHg19.over.chain \
  path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted_hg19.bed \
  path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted_unmapped.bed

echo "Original regions:" $(wc -l < path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted.bed)
echo "Lifted regions:" $(wc -l < path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted_hg19.bed)
echo "Unmapped regions:" $(wc -l < path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted_unmapped.bed)

# Extend TSS by ±1Mb
awk 'BEGIN { OFS="\t" } { $2 = $2 - 1000000; $3 = $3 + 1000000; print }' resources/Pilot_TSS_hg19.txt > path/to/output/directory/Pilot_TSS_hg19_1Mb_extended.bed

# Intersect hg19 ctcf sites with extended TSS cooridnates to extract all CTCF sites in the target gene loci
bedtools intersect \
    -a <(bedtools sort -i path/to/output/directory/ENCFF072BUT_HCT116_CTCF_sorted_hg19.bed | bedtools merge -i -) \
    -b <(bedtools sort -i path/to/output/directory/Pilot_TSS_hg19_1Mb_extended.bed | bedtools merge -i -) \
    -u > path/to/output/directory/CTCF_within_target_gene_loci_hg19.bed

# Make 4-column CTCF bed for bigWigAverageOverBed
awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $1":"$2"-"$3}' path/to/output/directory/CTCF_within_target_gene_loci_hg19.bed > path/to/output/directory/4col_CTCF_within_target_gene_loci_hg19.bed

# Generate an enhancers hg19 bed file
# Skip header with 'NR>1', split first column 'name' by ":" and "-", then print chr, start, end, name
awk -F',' 'NR>1 {
# Split the "name" column (first column) by "[:-]"
  split($1, a, "[:-]");
  print a[1], a[2], a[3], $1
}' OFS='\t' resources/CRUDO_FF_enhancers.csv > path/to/output/directory/hg19_enhancers_noTSS.bed


# Compute average BigWig signal over enhancers and CTCF sites
for bw in path/to/GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bw path/to/GSM2809610_Rao-2017-CHIP002-RAD21-treated.bw; do
    n=$(basename "$bw" | sed 's/.bw//')
    echo "Processing $n ..."
    path/to/bigWigAverageOverBed "$bw" path/to/output/directory/hg19_enhancers_noTSS.bed path/to/output/directory/CRUDO_hg19_enhancers_noTSS_signal_"$n".bed
    path/to/bigWigAverageOverBed "$bw" path/to/output/directory/4col_CTCF_within_target_gene_loci_hg19.bed path/to/output/directory/CRUDO_CTCF_within_target_gene_loci_hg19_signal_"$n".bed
done