import scanpy as sc
import stlearn as st
import matplotlib.pyplot as plt
import numpy as np
import os


plt.rcParams['figure.dpi'] = 200
plt.rcParams['savefig.dpi'] = 600
plt.rcParams['path.simplify'] = False

PBBCHF_DX = sc.read_h5ad("path/to/Sig_Simple_Grid_PBBCHF_DX.h5ad")
PBBCHF_PT = sc.read_h5ad( "path/to/Sig_Simple_Grid_PBBCHF_PT.h5ad")

best_lrs = ["CCL3_CCR5","CCL4_CCR3","CCL3_CCR4","CCL4_CCR5"]

out_dir = "path/to/LR_figures_dualpanel_obsm"
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
    vmax = np.percentile(combined, 99)  # robust scale


    fig, axes = plt.subplots(
        ncols=2,
        figsize=(16,8),
        gridspec_kw={'wspace':0.05}
    )

    # DX
    st.pl.lr_result_plot(
        PBBCHF_DX,
        use_result='lr_scores',
        use_lr=lr,
        show_color_bar=False,
        ax=axes[0],
        size=20,
        vmin=vmin,
        vmax=vmax
    )

    # PT
    st.pl.lr_result_plot(
        PBBCHF_PT,
        use_result='lr_scores',
        use_lr=lr,
        show_color_bar=False,
        ax=axes[1],
        size=20,
        vmin=vmin,
        vmax=vmax
    )


    for ax in axes:
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

    cbar = fig.colorbar(
        sm,
        ax=axes,
        fraction=0.025,
        pad=0.02
    )

    cbar.ax.tick_params(labelsize=12)

    plt.tight_layout()

    plt.savefig(
        f"{out_dir}/{lr}_DX_vs_PT_lr_scores.png",
        dpi=600,
        bbox_inches="tight"
    )

    plt.close(fig)

print("All LR-specific shared-scale dual-panel figures saved.")

