#!/bin/bash
#SBATCH --job-name=Simple_CCI_LOOP
#SBATCH --partition=ompi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --mem=20G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-199
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
set -x

####################################
# User settings
####################################
H5AD_DIR="/path/to/processed_data"
OUT_DIR="/path/to/stlearn/celltype/grid125/LR"
USE_LABEL="celltype"
GRID_N=125

mkdir -p logs

####################################
# Collect input files
####################################
mapfile -t H5AD_FILES < <(ls ${H5AD_DIR}/*_Anno.h5ad | sort)

N_FILES=${#H5AD_FILES[@]}

if [[ ${SLURM_ARRAY_TASK_ID} -ge ${N_FILES} ]]; then
    echo "Task ${SLURM_ARRAY_TASK_ID} >= ${N_FILES}, exiting."
    exit 0
fi

H5AD_FILE="${H5AD_FILES[$SLURM_ARRAY_TASK_ID]}"
SAMPLE_ID=$(basename "${H5AD_FILE}" _Anno.h5ad)

echo "Running sample: ${SAMPLE_ID}"
echo "Input file:   ${H5AD_FILE}"

####################################
# Environment
####################################
source ~/miniconda3/etc/profile.d/conda.sh
conda activate stlearn

export NUMBA_DISABLE_CACHE=1
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_PER_TASK

####################################
# Run Python
####################################
python << EOF
import os
import stlearn as st
import scanpy as sc

print("Sample:", "${SAMPLE_ID}")

adata = sc.read_h5ad("${H5AD_FILE}")

# -----------------------------
# QC
# -----------------------------
st.pp.filter_genes(adata, min_counts=10)
st.pp.filter_cells(adata, min_counts=10)

adata.X = adata.X.toarray()

# -----------------------------
# PCA + neighbors
# -----------------------------
st.em.run_pca(adata, n_comps=50, random_state=0)
st.pp.neighbors(adata, n_neighbors=25, use_rep="X_pca", random_state=0)
st.pp.normalize_total(adata)

# -----------------------------
# Grid
# -----------------------------
grid = st.tl.cci.grid(
    adata,
    n_row=${GRID_N},
    n_col=${GRID_N},
    use_label="${USE_LABEL}"
)

# -----------------------------
# Load LR database
# -----------------------------
lrs = st.tl.cci.load_lrs(
    ["connectomeDB2020_lit"],
    species="human"
)

bad_genes = grid.var_names[grid.var_names.str.contains("_")]
grid = grid[:, ~grid.var_names.isin(bad_genes)].copy()

# -----------------------------
# CCI rank
# -----------------------------
st.tl.cci.run(
    grid,
    lrs,
    min_spots=20,
    distance=250,
    n_pairs=10000,
    n_cpus=None
)

st.tl.cci.adj_pvals(
    grid,
    correct_axis="spot",
    pval_adj_cutoff=0.05,
    adj_method="fdr_bh"
)

# -----------------------------
# Save
# -----------------------------
out = "${OUT_DIR}/Simple_Grid_${SAMPLE_ID}.h5ad"
grid.write(out)

print("Saved:", out)
EOF

echo "Finished sample: ${SAMPLE_ID}"
