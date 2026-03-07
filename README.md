# Spatially resolved T cell receptor tracking links tumor-immune segregation to clonotype-specific T cell states in high-risk neuroblastoma

## Abstract
High-risk neuroblastoma (HRNB) is a leading cause of pediatric cancer death, and the contribution of T cells to disease biology and treatment response remains poorly defined. To address this, we combined longitudinal bulk T cell receptor (TCR) sequencing with spatial transcriptomics and in situ TCR profiling in paired pre- and post-therapy HRNB samples to map T cell geography and phenotype. We found that tumors were not uniformly immune-cold but instead contained distinct spatial immune niches and tumor-dense regions. Across therapy, samples diverged along two ecological trajectories yet consistently retained this tumor-immune segregation. Within that architecture, γδT cells accumulated in neuroblast-rich regions, and identical TCR clonotypes exhibited less differentiated phenotypes in immune-rich areas but dysfunction-associated programs at the tumor interface. Spatial context therefore emerges as a major determinant of T cell state in HRNB, nominating tumor-interface γδT cells and niche remodeling as potential therapeutic levers.

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

8. **cci/**  
   Cell-cell interaction analysis using stlearn.

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
