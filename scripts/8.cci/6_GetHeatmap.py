import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Try cross-platform font paths; fall back to default if Helvetica not found
font_paths = [
    (r"C:\Users\19808\AppData\Local\Microsoft\Windows\Fonts\Helvetica.ttf", "Helvetica"),
    ("/System/Library/Fonts/Helvetica.ttc", "Helvetica"),
    ("/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf", "Liberation Sans"),
]
font_added = False
for fp, family in font_paths:
    if os.path.exists(fp):
        fm.fontManager.addfont(fp)
        plt.rcParams["font.family"] = family
        font_added = True
        break
if not font_added:
    plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.size"] = 7


df = pd.read_csv(
    "LR_logFC_PT_vs_DX.tsv",
    sep="\t"
)


heatmap_df = df.pivot(
    index="LR",
    columns="Patient",
    values="logFC_PT_vs_DX"
)


fig_height = len(heatmap_df) * 0.15
fig, ax = plt.subplots(figsize=(5, fig_height))

hm = sns.heatmap(
    heatmap_df,
    cmap="bwr",
    center=0,
    linewidths=0.15,
    linecolor="white",
    cbar_kws={
        "label": "log2(PT/DX)",
        "shrink": 0.45,
        "aspect": 35,
        "pad": 0.015
    },
    ax=ax
)


ax.set_xlabel("Patient", fontsize=7)
ax.set_ylabel("Ligand-Receptor", fontsize=7)
ax.set_title("LR Network log2 Fold Change (PT vs DX)", fontsize=7)

ax.tick_params(axis="x", labelsize=7, rotation=45)
ax.tick_params(axis="y", labelsize=7)


cbar = hm.collections[0].colorbar
cbar.ax.tick_params(labelsize=7)
cbar.set_label("log2(PT/DX)", fontsize=7)


plt.tight_layout()

plt.savefig(
    "LR_network_logFC_PT_vs_DX_heatmap.pdf",
    format="pdf"
)

plt.close()