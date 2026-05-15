# Spatially resolved T cell receptor tracking links tumor-immune segregation to clonotype-specific T cell states in high-risk neuroblastoma

## Abstract
High-risk neuroblastoma (HRNB) is a leading cause of pediatric cancer death. Current therapies center on intensive multimodal treatment including anti-GD2 therapy, with growing interest in harnessing T cell-mediated immunity. How T cells and their receptors (T-cell receptors, TCRs) are spatially organized and function within tumors remains poorly defined. To assess whether intratumoral location influences clonotype-specific T cell states, we profiled TCR repertoires across blood and tumor samples from 37 patients with HRNB using longitudinal bulk TCR sequencing. In a nested subset of 5 patients with paired pre- and post-therapy tumors, we integrated spatial transcriptomics with in situ TCR profiling. Across all tumors, T and B cells preferentially co-localized in immune-rich regions and showed reduced proximity to neuroblast cells. Despite this compartmentalized architecture, γδT cells were more evenly distributed across tumor sections and showed greater proximity to neuroblast-rich regions than other T cell subsets. Within TCR clonotypes, spatial location was associated with distinct transcriptional states, with immune-rich regions supporting more progenitor-like programs. These findings identify spatial context as a key determinant of phenotype clonotype-specific T cell phenotype and highlight γδT cells cells as a spatially distinct population with potential roles in neuroblastoma tumor-immune interactions.

## Repository Structure

```
neuroblastoma_spatial_tcr/
├── README.md
├── requirements.txt
├── manuscript_sourcedata/ # source data in Excel format
├── data/
│   ├── README.md
│   └── processed/ # source data for figure generation in .csv format
└── scripts/
    ├── 1.process_and_integrate.ipynb
    ├── 2.subcluster_and_visualize_metadata.ipynb
    ├── 3.calculate_celltype_proportions_and_neighborhoods.ipynb
    ├── 4.correlate_distance_and_genelists.ipynb
    ├── 5.find_domains.ipynb
    ├── 6.tcr_bulk/
    │   ├── 1_run_mixcr_command_example.sh
    │   ├── 2_immunarch_tcrs.R
    │   └── 3_make_sankey_plot.ipynb
    └── 7.find_tcr_clonotypes.ipynb
```

## Overview of Notebooks and Scripts

1. **process_and_integrate**  
   Preprocessing, integration, and harmonization of spatial transcriptomics 
data.

2. **subcluster_and_visualize_metadata**  
   Subclustering and metadata-based visualization of spatially resolved cell 
populations.

3. **calculate_celltype_proportions_and_neighborhoods**  
   Computes cell type proportions and spatial neighborhood 
metrics/enrichment.

4. **correlate_distance_and_genelists**  
   Spatial distance calculations and correlation with gene signatures.

5. **find_domains**  
   Identification of latent spatial domains using scimap.

6. **tcr_bulk/**  
   CapTCR-seq processing and clonotype anlaysis pipeline including MiXCR4, 
Immunarch and TCR probe visualization.

7. **find_tcr_clonotypes**  
   Clonotype-aware spatial analysis.

## Software Requirements

### Python Environment
See `requirements.txt` for Python package dependencies.

### External Tools

**MiXCR** v4.6.0 (https://github.com/milaboratory/mixcr)

**Xenium Explorer** v4.1

**R** v4.4.0
- immunarch v0.10.3
- ggplot2 v4.0.0
- dplyr v1.1.4

## Reproducibility

This repository supports two levels of reproducibility.

1. **Full analysis reproduction from the original `.h5ad` file.**  
   The `.h5ad` file is not included because it contains controlled-access patient-derived data that will be used in a future manuscript, but it is available from the authors upon reasonable request.

2. **Figure reproduction from processed CSV files.**  
   All processed `.csv` files required to regenerate the manuscript figures are provided in `data/processed/`.

The manuscript source data Excel files are provided in `manuscript_sourcedata/`. These files are intended for journal source-data submission and reader inspection. The canonical code-readable files are the `.csv` files in `data/processed/`.
