# Neuroblastoma Spatial TCR Analysis

This repository contains custom code for Xenium and CapTCR-seq TCR profiling in neuroblastoma samples.

## Repository Structure

```
neuroblastoma_spatial_tcr/
├── scripts/
│   ├── 1.process_and_integrate.ipynb
│   ├── 2.subcluster_and_visualize_metadata.ipynb
│   ├── 3.calculate_celltype_proportions_and_neighborhoods.ipynb
│   ├── 4.correlate_distance_and_genelists.ipynb
│   ├── 5.find_domains.ipynb
│   ├── 6.tcr_bulk/
│   │   ├── 1_run_mixcr_command_example.sh
│   │   ├── 2_immunarch_tcrs.R
│   │   └── 3_make_sankey_plot.ipynb
│   ├── 7.find_tcr_clonotypes.ipynb
│   └── 8.cci/
│       ├── 1_AddCellType_toRawObj.py
│       ├── 2_Run_stlearn_LR_Analysis.slurm
│       ├── 3_Run_CCI_Analysis.slurm
│       ├── 4_GetLRSummaryPerSample.py
│       ├── 5_GetSpecificLRContribution.py
│       ├── 6_GetHeatmap.py
│       ├── 7_GetSpatialMap_IndependentColorBar.py
│       ├── 8_Barplots.py
│       └── 9_Dotplot_BySingleCell.py
│
├── manu_figs/            # manuscript figures
│
└── manu_sourcedata/      # source data tables associated with manuscript figures
```

## Notebook Overview

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

8. **cci**
   Cell-cell interaction analysis using stlearn.

## Software Requirements

### Python Environment
See `requirements.txt` for Python package dependencies.

### External Tools

**MiXCR** v4.6.0 (https://github.com/milaboratory/mixcr)

**R** v4.4.0
- immunarch v0.10.3
- ggplot2 v4.0.0
- dplyr v1.1.4
