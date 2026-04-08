# Spatially resolved T cell receptor tracking links tumor-immune segregation to clonotype-specific T cell states in high-risk neuroblastoma

## Abstract
High-risk neuroblastoma (HRNB) is a leading cause of pediatric cancer death. Current therapies center on intensive multimodal treatment with anti-GD2 therapy, with growing interest in harnessing T cell-mediated immunity. However, how T cells and their receptors (T-cell receptors; TCRs) are spatially organized and function within tumors remains poorly defined. To assess whether intratumoral location influences clonotype-specific T cell states, we combined longitudinal bulk TCR sequencing in blood and tumor with spatial transcriptomics and in situ TCR profiling in paired pre- and post-therapy HRNB samples. Across all tumors, spatial compartmentalization is consistently present, with T and B cells segregated from neuroblasts. Despite this segregated architecture, γδT cells preferentially accumulated in neuroblast-rich regions, positioning them as the dominant T cell subset at the tumor interface. Within clonotypes, spatial location was associated with distinct transcriptional states, with immune-rich regions supporting more progenitor-like programs. These findings identify spatial context as a key determinant of phenotype clonotype-specific T cells, and nominate γδT cells as potential therapeutic levers.

## Repository Structure

```
neuroblastoma_spatial_tcr/
├── README.md
├── requirements.txt
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
