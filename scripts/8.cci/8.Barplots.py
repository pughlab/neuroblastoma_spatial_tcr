# Explore the cell-type composition of all spots and hotspots with high ligand–receptor scores.
# Mean cell-type proportions across the five samples were calculated to represent each group.
"""
1. For each LR interaction: Extract LR interaction scores from each spatial grid spot; Identify hotspot spots defined as the top p proportion (default: top 10%) of LR scores.

2. For each sample: Retrieve cell-type mixture proportions per spot; Convert proportions to weighted cell counts using the number of cells per spot; Calculate cell-type composition for all spots (ALL) and hotspot spots (HOT)

3. For each condition (DX or PT): Aggregate results across all samples; Compute the mean cell-type composition for: DX_all, DX_hot, PT_all,PT_hot

4. Visualization: For each LR interaction, generate a stacked bar plot showing the mean weighted cell-type composition
     for the four conditions above.

"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['path.simplify'] = False


lrs = ["CCL4_CCR5","CCL4_CCR3","CCL3_CCR4","CCL3_CCR5"]
p = 0.1

samples = {
"PBBCHF_DX": sc.read_h5ad("path/to/Sig_Simple_Grid_PBBCHF_DX.h5ad"),
"PBBCHF_PT": sc.read_h5ad("path/to/Sig_Simple_Grid_PBBCHF_PT.h5ad"),

"PAZNRG_DX": sc.read_h5ad("path/to/Sig_Simple_Grid_PAZNRG_DX.h5ad"),
"PAZNRG_PT": sc.read_h5ad("path/to/Sig_Simple_Grid_PAZNRG_PT.h5ad"),

"PAZWZN_DX": sc.read_h5ad("path/to/Sig_Simple_Grid_PAZWZN_DX.h5ad"),
"PAZWZN_PT": sc.read_h5ad("path/to/Sig_Simple_Grid_PAZWZN_PT.h5ad"),

"PBADJC_DX": sc.read_h5ad("path/to/Sig_Simple_Grid_PBADJC_DX.h5ad"),
"PBADJC_PT": sc.read_h5ad("path/to/Sig_Simple_Grid_PBADJC_PT.h5ad"),

"PBBKFP_DX": sc.read_h5ad("path/to/Sig_Simple_Grid_PBBKFP_DX.h5ad"),
"PBBKFP_PT": sc.read_h5ad("path/to/Sig_Simple_Grid_PBBKFP_PT.h5ad")
}

# -------------------------
# celltype color
# -------------------------
color_map = {
'Endothelial': '#d73027',
'Fibroblast': '#f46d43',
'Schwann': '#fdae61',
'Neuroblast': '#fee090',
'Macrophage': '#e0f3f8',
'B': '#abd9e9',
'T': '#74add1'
}

out_dir = "LR_group_mean_barplots"
os.makedirs(out_dir, exist_ok=True)


for lr in lrs:

    dx_all_list = []
    dx_hot_list = []
    pt_all_list = []
    pt_hot_list = []

    for name, adata in samples.items():
        lr_matches = np.where(adata.uns["lr_summary"].index.values == lr)[0]
        if len(lr_matches) == 0:
            print(f"Skipping {lr}: not found in lr_summary for {name}")
            continue
        lr_index = lr_matches[0]

        vals = adata.obsm["lr_scores"][:, lr_index]

        cutoff = np.quantile(vals, 1 - p)
        mask = vals >= cutoff

        cell_counts = adata.obs["n_cells"]
        mix_prop = adata.uns["celltype"]

        weighted = mix_prop.multiply(cell_counts, axis=0)

        mix_all_sum = weighted.sum()
        mix_all_prop = mix_all_sum / mix_all_sum.sum()

        mix_hot_sum = weighted.loc[mask].sum()
        mix_hot_prop = mix_hot_sum / mix_hot_sum.sum()

        if "_DX" in name:
            dx_all_list.append(mix_all_prop)
            dx_hot_list.append(mix_hot_prop)

        if "_PT" in name:
            pt_all_list.append(mix_all_prop)
            pt_hot_list.append(mix_hot_prop)

    if not dx_all_list or not pt_all_list:
        print(f"Skipping {lr}: no samples with this LR in lr_summary")
        continue

    dx_all_mean = pd.DataFrame(dx_all_list).fillna(0).mean()
    dx_hot_mean = pd.DataFrame(dx_hot_list).fillna(0).mean()
    pt_all_mean = pd.DataFrame(pt_all_list).fillna(0).mean()
    pt_hot_mean = pd.DataFrame(pt_hot_list).fillna(0).mean()

    plot_df = pd.DataFrame({
        "DX_all": dx_all_mean,
        "DX_hot": dx_hot_mean,
        "PT_all": pt_all_mean,
        "PT_hot": pt_hot_mean
    }).fillna(0)

    plot_df = plot_df.sort_index()

    ordered_cols = [c for c in color_map if c in plot_df.index]
    plot_df = plot_df.loc[ordered_cols]
    colors = [color_map[c] for c in ordered_cols]


    plt.figure(figsize=(8,6))

    plot_df.T.plot(
        kind="bar",
        stacked=True,
        color=colors,
        edgecolor="none"
    )

    plt.ylabel("Mean weighted cell composition")
    plt.title(f"{lr} composition (Top {int(p*100)}% hotspot)")
    plt.xticks(rotation=0)
    plt.legend(bbox_to_anchor=(1.05,1))
    plt.tight_layout()

    plt.savefig(
        f"{out_dir}/{lr}_DX_PT_mean_barplot.pdf"
    )

    plt.close()

print("All LR mean barplots saved.")