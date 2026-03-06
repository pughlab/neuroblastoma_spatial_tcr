import scanpy as sc
import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import os


plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['path.simplify'] = False

PBBCHF_DX = sc.read_h5ad("path/to/Sig_Simple_Grid_PBBCHF_DX.h5ad")
PBBCHF_PT = sc.read_h5ad("path/to/Sig_Simple_Grid_PBBCHF_PT.h5ad")

best_lrs = ["CCL3_CCR5","CCL4_CCR3","CCL3_CCR4","CCL4_CCR5"]

out_dir = "path/to/LR_figures_singlepanel"
os.makedirs(out_dir, exist_ok=True)


for lr in best_lrs:
    lr_matches = np.where(PBBCHF_DX.uns["lr_summary"].index.values == lr)[0]
    if len(lr_matches) == 0:
        print(f"Skipping {lr}: not found in lr_summary")
        continue
    lr_index = lr_matches[0]

    dx_vals = PBBCHF_DX.obsm["lr_scores"][:, lr_index]
    pt_vals = PBBCHF_PT.obsm["lr_scores"][:, lr_index]

    combined = np.concatenate([dx_vals, pt_vals])

    vmin = 0
    vmax = np.percentile(combined, 99)

    for adata, suffix in [(PBBCHF_DX, "PBBCHF_DX"), (PBBCHF_PT, "PBBCHF_PT")]:
        fig, ax = plt.subplots(figsize=(8, 8))

        st.pl.lr_result_plot(
            adata,
            use_result='lr_scores',
            use_lr=lr,
            show_color_bar=False,
            ax=ax,
            size=20,
            vmin=vmin,
            vmax=vmax
        )

        for coll in ax.collections:
            coll.set_edgecolor("none")
            coll.set_linewidth(0)
            coll.set_antialiased(False)

        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel("")

        sm = plt.cm.ScalarMappable(
            cmap=plt.cm.Spectral_r,
            norm=plt.Normalize(vmin=vmin, vmax=vmax)
        )
        sm.set_array([])

        cbar = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.ax.tick_params(labelsize=12)

        plt.tight_layout()
        plt.savefig(f"{out_dir}/{lr}_{suffix}_lr_scores.png", dpi=600, bbox_inches="tight")
        plt.close(fig)

print("All LR figures saved.")