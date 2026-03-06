#5. Calculate contribution of specific LR interaction to sample level CCI network and log2 fold change between PT and DX
"""
This script processes the merged ligand–receptor (LR) interaction summary table
generated from stLearn CCI analysis and computes normalized LR interaction
proportions as well as log2 fold changes between time points.

Steps performed in this analysis:

1. Load the merged LR summary table containing LR interaction statistics for
   multiple samples. For each sample, calculate the total number of significant cell-type-level
   LR interaction spots:
        total_CCI = sum(n-spot_cci_sig_celltype)

2. Normalize LR interaction counts within each sample by computing the network
   proportion of each LR pair:
       LR_prop_network = n-spot_cci_sig_celltype / total_CCI

   This represents the fraction of all LR interactions within a sample that are attributed to a specific LR pair.

4. Compute the log2 fold change in LR interaction proportion between PT and DX:
       logFC_PT_vs_DX = log2((PT + ε) / (DX + ε))

   A small constant (ε = 1e-10) is added to avoid division by zero.
"""

import pandas as pd
import numpy as np

lr = pd.read_csv(
    "path/to/All_LR_Summary.tsv",
    sep="\t"
)


print(lr["LR"].unique())

lr_df = lr.copy()

lr_df["total_CCI"] = (
    lr_df.groupby("Sample")["n-spot_cci_sig_celltype"]
    .transform("sum")
)

lr_df["LR_prop_network"] = (
    lr_df["n-spot_cci_sig_celltype"] /
    lr_df["total_CCI"]
)

logFC_df = (
    lr_df[["LR", "Patient", "TimePoint", "LR_prop_network"]]
    .pivot_table(
        index=["LR", "Patient"],
        columns="TimePoint",
        values="LR_prop_network"
    )
    .reset_index()
)

# Handle missing PT or DX columns (e.g. when a patient has only one timepoint)
eps = 1e-10
pt_vals = logFC_df["PT"].fillna(eps).values if "PT" in logFC_df.columns else np.full(len(logFC_df), eps)
dx_vals = logFC_df["DX"].fillna(eps).values if "DX" in logFC_df.columns else np.full(len(logFC_df), eps)
logFC_df["logFC_PT_vs_DX"] = np.log2((pt_vals + eps) / (dx_vals + eps))


print(logFC_df.head())

logFC_df.to_csv(
    "LR_logFC_PT_vs_DX.tsv",
    sep="\t",
    index=False
)