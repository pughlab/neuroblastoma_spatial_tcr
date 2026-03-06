#

#### color map
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.font_manager as fm

font_path = r"path\to\Helvetica.ttf"
fm.fontManager.addfont(font_path)

plt.rcParams["font.family"] = "Helvetica"
#plt.rcParams["font.size"] =7
#plt.rcParams["axes.titlesize"] = 7
#plt.rcParams["axes.labelsize"] = 7
#plt.rcParams["xtick.labelsize"] =7
#plt.rcParams["ytick.labelsize"] = 7





adata=sc.read_h5ad("path/to/xenium_integrated_labeled.h5ad")



color_map = {
    'Endothelial': '#d73027',
    'Fibroblast': '#f46d43',
    'Schwann': '#fdae61',
    'Neuroblast': '#fee090',
    'Macrophage': '#e0f3f8',
    'B': '#abd9e9',
    'T': '#74add1',
}


adata.obs["celltype"] = adata.obs["celltype"].astype("category")

ordered = [c for c in color_map if c in adata.obs["celltype"].cat.categories]
adata.obs["celltype"] = adata.obs["celltype"].cat.reorder_categories(ordered)

adata.obs["timepoint"].unique()

adata_dx = adata[adata.obs["timepoint"] == "DX"].copy()
adata_pt = adata[adata.obs["timepoint"] == "PT"].copy()

genes = ["CCL3","CCL4","CCR3","CCR4","CCR5"]


# DX
dp = sc.pl.dotplot(
    adata_dx,
    var_names=genes,
    groupby="celltype",
    swap_axes=True,
    standard_scale="var",
    cmap="YlGnBu",
    return_fig=True,
    show=False
)

dp.savefig("dotplot_DX.pdf")


# PT
dp = sc.pl.dotplot(
    adata_pt,
    var_names=genes,
    groupby="celltype",
    swap_axes=True,
    standard_scale="var",
    cmap="YlGnBu",
    return_fig=True,
    show=False
)

dp.savefig("dotplot_PT.pdf")