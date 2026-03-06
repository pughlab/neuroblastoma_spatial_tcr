#1. Add Cell Type Annotation to Raw Adata Object
"""
st.tl.cci.grid in stLearn requires adata.uns["spatial"] (image metadata), which is missing in the processed merged AnnData object containing the 10 samples.
Therefore, each sample AnnData was reconstructed from the raw Xenium output directories and the necessary obs annotations were re-added.
"""

###Create Obj With st.ReadXenium
#Use PBBKFP_DX output directory as example


import os
os.environ["NUMBA_DISABLE_CACHE"] = "1"
import stlearn as st
import scanpy as sc

adata = st.ReadXenium(feature_cell_matrix_file= "cell_feature_matrix.h5",
                      cell_summary_file= "cells.csv.gz",
                      library_id= "PBBKFP_PT", #
                      #image_path= "morphology.ome.tif",
                      #scale=1,
                      spot_diameter_fullres=15,
                      #alignment_matrix_file=data_dir / "he_imagealignment.csv",
                      experiment_xenium_file="experiment.xenium" )


adata.write("./PBBKFP_PT.h5ad")



###Add cell type annotation to raw object from ReadXenium

raw = sc.read_h5ad("path/to/Read Xenium created object/PBBKFP_PT.h5ad")
anno = sc.read_h5ad("path/to/fully annotated object/xenium_integrated_labeled.h5ad")

#Subset corresponding sample from integrated h5ad
anno = anno[anno.obs["patient"] == "PBBKFP"].copy()
anno_sub = anno[anno.obs["sampleName"] == "PBBKFP_PT"].copy()
anno_obs_by_cellid = anno_sub.obs.copy()
anno_obs_by_cellid.index = anno_obs_by_cellid["cell_id"]

common_cells = raw.obs_names.intersection(anno_obs_by_cellid.index)
print(len(common_cells), raw.n_obs)

raw_filt = raw[common_cells].copy()
print(len(common_cells), raw_filt.n_obs)


raw_filt.obs = raw_filt.obs.join(anno_obs_by_cellid)


raw_filt.write("./PBBKFP_PT_Anno.h5ad")

