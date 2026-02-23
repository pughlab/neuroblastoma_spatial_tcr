# Neuroblastoma Spatial TCR Analysis

This repository contains custom code for Xenium and CapTCR-seq TCR profiling in neuroblastoma samples.

## Repository Structure

scripts/
- 1.process_and_integrate.ipynb  
- 2.subcluster_and_visualize_metadata.ipynb  
- 3.calculate_celltype_proportions_and_neighborhoods.ipynb
- 4.correlate_distance_and_genelists.ipynb  
- 5.find_domains.ipynb  

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
   Identification of latent spatial domains.


## Bulk TCR Analysis Workflow

Directory: `scripts/6.tcr_bulk/`

This directory contains the bulk TCR processing and clonotype analysis 
pipeline:

1. **MiXCR**
   - `1_run_mixcr_command_example.sh`
   - Alignment and clonotype assembly from raw TCR sequencing data.

2. **GLIPH2**
   - `2_concatenate_gliph_input.R`
   - `3_gliph_community_detection.R`
   - Clonotype clustering and antigen-specific community detection.

3. **Immunarch**
   - `4_immunarch_tcrs.R`
   - Diversity, repertoire metrics, and comparative analysis.

4. **TCR probe visualization**
   - `5_make_sankey_plot.ipynb`
   - Sankey-based visualization of clonotype relationships and dynamics.

Pipeline summary:

MiXCR → GLIPH2 → Immunarch → Visualization

## Environment

Coming soon.

