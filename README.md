# CRUDO analysis workflow
This repository provides code for analyzing CRISPRi of regulatory elements upon degron operation (CRUDO) datasets — a method for high-throughput enhancer perturbation with acute, auxin-inducible depletion of a protein of interest (in this case, cohesin).

For more details, see our manuscript [3D contacts drive enhancer-promoter regulation through cohesin-dependent and -independent mechanisms](https://www.biorxiv.org/content/10.1101/2024.07.12.603288v2). More experimental details can be found on our [protocols.io](https://www.protocols.io/view/crispri-of-regulatory-elements-upon-degron-operati-dfbk3ikw.html). 

## Requirements
### Software requirements
  - macOS 14.7.3 (tested)

### Python requirements
  -  Python version: 3.9.23 (tested)
  -  Dependencies are listed in requirements.txt.

## Data availability
All datasets used in the manuscript are listed in the Usage section.
Some are included in the `resources/` folder.
If using alternative datasets, update script paths and parameters accordingly.

## Preprocessing tools
Data were preprocessed with:
  -  https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction
  -  https://github.com/karbalayghareh/ENCODE-rE2G
  -  https://github.com/EngreitzLab/crispri-flowfish
  -  https://github.com/argschwind/TAPseq
  -  https://github.com/argschwind/TAPseq_workflow
  -  https://github.com/ENCODE-DCC/chip-seq-pipeline2 (H3K27ac BAM files)

## Setting up the environment
Follow these steps to set up an environment with all required Python packages:

1. Create a virtual environment
``` bash
python3.9 -m venv venv
```
2. Activate the environment
``` bash
source venv/bin/activate
```
3. Install the dependencies
``` bash
pip install -r PackageRequirements.txt
```

## Usage
Follow the steps below to run the CRUDO analysis pipeline.
Each script has detailed docstrings describing its inputs and outputs.

### Step 1: Extract Hi-C contacts
Data sources: 
- Hi-C files can be downloaded from ENCODE:
    - untreated: https://www.encodeproject.org/files/ENCFF528XGK/
    - treated: https://www.encodeproject.org/files/ENCFF317OIA/
- The CRUDO target list is provided in `resources/`.

Run:
``` bash
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
```

Note: If running multiple times with different signal types/normalizations, merge the resulting CSVs separately.

### Step 2: Compute element features
Provided in `resources/`:
  - `TargetList.csv` (CRUDO target list)
  - `df_ABC_pilot.csv` (unthresholded ABC scores for pilot genes untreated and treated)
  - `ENCFF072BUT_HCT116_CTCF.bed` (CTCF peaks in the cohesin presence condition)

H3K27ac BAM files can be generated using the ENCODE ChIP-seq pipeline:
  - H3K27ac untreated: SRX3275590 (input: SRX4336996)
  - H3K27ac treated: SRX3275591 (input: SRX4336997)
(alternatively contact author for files too large to host)

Run: 
``` bash
    python scripts/_02_get_element_features.py \
        --element_file resources/TargetList.csv \
        --abc_file resources/df_ABC_pilot.csv \
        --ctcf_bed resources/ENCFF072BUT_HCT116_CTCF.bed  \
        --h3k27ac_untreated_bam path/to/SRR6164278_H3K27Ac-untreated.srt.nodup.bam \
        --h3k27ac_treated_bam path/to/SRR6164279_H3K27Ac-treated.srt.nodup.bam \
        --output_directory path/to/output/directory/
```

### Step 3: Compute CRUDO element effect sizes and cohesin dependency for FlowFISH data
Provided in `resources/`:
  - `TargetList.csv` (CRUDO target list)
  - `byExperimentRep/` (processed guide-level CRUDO-FF screen data)
    
Default: --alpha_differential 0.05
    (In our analysis, this was adjusted to 0.0020833333333333333, derived from running `04b_CRUDO_FF_enhancer_classification.py`)

Run:
``` bash
    python scripts/_03_CRUDO_FF_analysis_pipeline.py  \
        --element_file resources/TargetList.csv \
        --processed_FF_directory resources/byExperimentRep/ \
        --output_directory path/to/output/directory/ \
        --alpha_differential 0.0020833333333333333
```

### Step 4: Enhancer classification
#### 1. Merge output CSVs from Steps 1–3
   
Run:
``` bash
    python scripts/_04a_CRUDO_FF_datasets_merge.py \
        --hic_path path/to/TargetList.SCALE.normalized.merged.5Kb.csv\
        --features_path path/to/output/directory/CRUDO_element_features.csv \
        --elements_effects_path path/to/output/directory/CRUDO_FF_ElementLevel_analysis.csv\
        --output_directory path/to/output/directory/
```

#### 2. Classify enhancers:

input = merged file from Step 4.1, which is the equivalent to Supplementary Table 2c in our manuscript and is provided in `resources/`

Run:
``` bash
    python scripts/_04b_CRUDO_FF_enhancer_classification.py \
        --input_file resources/CRUDO_FF_SupplementaryTable2c.csv \
        --output_directory path/to/output/directory/
```
### Step 5: Analyze RAD21 retention
#### 1. Prepare bed files

First, adjust file paths in the script to match your environment.

Required output directory (adjust as needed):
  - `path/to/output/directory/`

Provided in `resources/`:
  - `ENCFF072BUT_HCT116_CTCF.bed`
  - `Pilot_TSS_hg19.txt`
  - `CRUDO_FF_enhancers.csv` (ouput from Step #4)

To download from GEO (RAD21 signal files):
  - `path/to/GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bw`
  - `path/to/GSM2809609_Rao-2017-CHIP001-RAD21-untreated.bw`

To download from UCSC:
  - `path/to/liftOver`
  - `path/to/hg38ToHg19.over.chain`
  - `path/to/bigWigAverageOverBed`

Run:
``` bash
    bash scripts/_05a_PrepBeds_CRUDO_RAD21_analysis.sh
```

#### 2. Run the Analysis

The paths to the prepared BEDs are hardcoded in the analysis script. Update them to point to the files generated in Step 5.1.

Run:
``` bash
    python scripts/_05b_CRUDO_RAD21_analysis.py \
        --enhancers resources/CRUDO_FF_enhancers.csv \
        --output_directory_csv path/to/output/directory/ \
        --output_directory_plot path/to/output/directory/plots/
```

### Step 6: Visualize CRUDO enhancers

The input for this script is the output from step 5.2 `CRUDO_FF_enhancers_RAD21.csv` and is also provided in `resources/`

Run:
``` bash
    python scripts/_06_CRUDO_Enhancer_ElementLevel_Visualization.py \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --output_directory_plot path/to/output/directory/plots/
```

### Step 7: Model cohesin dependence for CRUDO enhancers
The input for this script is the output from step 5.2 `CRUDO_FF_enhancers_RAD21.csv` and is also provided in `resources/`

Run:
``` bash
    python scripts/_07_CRUDO_Enhancer_CohesinDep_modeling.py \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --output_directory_plot path/to/output/directory/plots/
```

### Step 8: Visualize guide-level effects at CRUDO enhancers

Provided in `resources/`:
  - `CRUDO_FF_enhancers_RAD21.csv` (CRUDO enhancers, ie output from step 5.2)
  - `byExperimentRep/` (processed guide-level CRUDO-FF screen data)

Run:
``` bash 
    python scripts/_08_CRUDO_FF_GuideLevel_Visualization.py  \
        --enhancers resources/CRUDO_FF_enhancers_RAD21.csv \
        --processed_FF_directory resources/byExperimentRep/ \
        --output_directory_plot path/to/output/directory/plots/
```      

### Step 9: Compute CRUDO element effect sizes and cohesin dependency for TAP-seq data

The processed guide-level CRUDO-DC-TAP-seq data is provided in `resources/`.

Run:
``` bash
    python scripts/_09_CRUDO_TAP_analysis_pipeline.py \
        --processed_guides_file resources/CRUDO_tapseq.results.guide_level.txt \
        --output_directory path/to/output/directory/
```     

### Step 10: Compare CRUDO element effects from FlowFISH and DC-TAP-seq

Provided in `resources`:
  - `CRUDO_FF_SupplementaryTable2c.csv` (CRUDO-FF element level effects, merged file from Step 4.1, which is the equivalent to Supplementary Table 2c in our manuscript)
  - `CRUDO_TAP_ElementLevel_analysis.csv` (CRUDO-DC-TAP-seq element level effects, output from Step 9)

Run:
``` bash
    python scripts/_10_CRUDO_TAPvFF_ElementLevel_Visualization.py \
        --FF_effects resources/CRUDO_FF_SupplementaryTable2c.csv \
        --TAP_effects resources/CRUDO_TAP_ElementLevel_analysis.csv \
        --output_directory_plot path/to/output/directory/plots/
```

### Step 11: Simulate the CRUDO experimental framework

This script runs with the following default parameters (can be changed if required):
  - `mRNA half-life in minutes: 600`
  - `Auxin removal start time in minutes: 30`
  - `dCas9 induction rate in minutes: 1440`
  - `Auxin induction time point in minutes: 1080`
  - `dCas9 induction time point in minutes: 0`
  - `Measurement time point in minutes: 1440`

The following parameters must be specified by the user. Example values:
  - `Cohesin effect on gene expression: 0.4`
  - `Enhancer effect on gene expression: 0.25`
  - `Cohesin effect on enhancer function: 0.9 `

Run:
``` bash
    python scripts/_11_CRUDO_FrameworkSimulation.py \
	    --ce 0.4 \
	    --e 0.25 \
	    --ae 0.9 \
	    --output_directory_plot path/to/output/directory/plots/
```

### Step 12: Compare CRUDO enhancer effects to published enhancer sets

Provided in `resources`:
  - `CRUDO_FF_SupplementaryTable2c.csv` (CRUDO-FF element level effects, merged file from Step 4.1, which is the equivalent to Supplementary Table 2c in our manuscript)
  - `Fulco2019_S3a.txt` (Fulco et al. 2019 Tiled CRISPRi screen enhancer effect from Supplementary Table 3a)
  - `Fulco2019_S6a.txt` (Fulco et al. 2019 Tiled CRISPRi screen TSS positions from Supplementary Table 6a)

To download from github (ENCODE rE2G benchmarking set of compiled CRISPRi screen results):
  - https://github.com/EngreitzLab/CRISPR_comparison/blob/main/resources/crispr_data/EPCrisprBenchmark_ensemble_data_GRCh38.tsv.gz


Note: Paths to published enhancer sets are hardcoded in the script and need to be adjusted accordingly.

Run:
``` bash
    python scripts/_12_CRUDOvPublishedEnhancers.py \
        --CRUDO_elements resources/CRUDO_FF_SupplementaryTable2c.csv \
        --output_directory_plot path/to/output/directory/plots/
```

### Step 13: Convert PRO-seq RPKM to TPM and visualize gene expression.

This script takes PRO-seq RPKM data from Rao et al., 2017 (`GSE106886_Rao-2017-Genes.rpkm.txt`), provided in `resources/`.

Run:
``` bash
    python scripts/_13_CRUDO_PROseq_GeneExpression.py \
        --PRO_RPKM resources/GSE106886_Rao-2017-Genes.rpkm.txt \
        --output_directory path/to/output/directory/
```

### Step 13: Analyze ABC enhancer-gene predictions with/without cohesin and stratify by CRUDO expression

Provided in `resources/`:
  - `1a_ABC_thresholded_noAux.csv` (thresholded ABC predictions in the presence of cohesin, which is the equivalent to Supplementary Table 1a in our manuscript)
  - `1b_ABC_thresholded_Aux.csv` (thresholded ABC predictions in the absence of cohesin, which is the equivalent to Supplementary Table 1b in our manuscript)
  - `CRUDO_Genes_Pro.tpm.txt` (PRO-seq TPM data (converted from Rao et al. 2017 RPKM data), output from step 13)
  - `Auxin_vs_Control.RAD21.Genes.DESeq2.txt` (Rao et al. 2017 DeSeq2 data)
  - `HousekeepingGenes_Fulco2019_S5d.txt` (A list of housekeeping genes defined in Fulco et al. 2019)
  - `CRUDO_FF_enhancers_RAD21.csv` (CRUDO enhancers, ie output from step 5.2)

Run:
``` bash
    python scripts/_14_ABC_analysis_CRUDO_stratification.py \
        --predictions_noAux  resources/1a_ABC_thresholded_noAux.csv\
        --predictions_Aux  resources/1b_ABC_thresholded_Aux.csv\
        --PRO_TPM resources/CRUDO_Genes_Pro.tpm.txt\
	    --PRO_DeSeq2 resources/Auxin_vs_Control.RAD21.Genes.DESeq2.txt\
        --housekeeping_genes  resources/HousekeepingGenes_Fulco2019_S5d.txt\
	    --CRUDO_enhancers  resources/CRUDO_FF_enhancers_RAD21.csv\
        --output_directory_plot path/to/output/directory/plots/
```

### Step 15: Plot Hi-C data
#### 1. Comvert .hic files to .cool

Hi-C files can be downloaded from ENCODE:
	- untreated: https://www.encodeproject.org/files/ENCFF528XGK/
 	- treated: https://www.encodeproject.org/files/ENCFF317OIA/

Run:
``` bash
    python scripts/_15a_Hic2Cool.py \
        --hic path/to/.hic \
        --output_directory path/to/output/directory/ \
        --resolution 5000 \
        --signal_type observed \
        --norm_method SCALE
```
#### 2. Plot enhancer-gene Hi-C squares

Use the .cool file generated in step 15.1.
Provided in `resources`:
	- `CRUDO_EnhancerGenePairs.csv` (genomic positions of CRUDO enhancer-gene pairs)

Run:
``` bash
    python _15b_CRUDO_PlotBinsHic.py \
        --cool path/to/.cool \
        --pairs CRUDO_EnhancerGenePairs.csv \
        --window-bp 20000 \
        --outdir path/to/output/directory/
```
