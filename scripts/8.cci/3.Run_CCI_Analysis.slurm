#!/bin/bash
#SBATCH --job-name=NoMix_CCI_run_cci_GRID125_CELLTYPE
#SBATCH --partition=ompi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=35
#SBATCH --mem=20G
#SBATCH --time=4-00:00:00

#SBATCH --array=0-13
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail
set -x

####################################
# User settings
####################################
GRID_DIR="/path/to/stlearn/celltype/grid125/LR"
OUT_DIR="/path/to/stlearn/celltype/grid125/NoMixture_CCI"

USE_LABEL="celltype"

CELL_PROP_CUTOFF=0.1
MIN_SPOTS=2
N_PERMS=1000

mkdir -p logs
mkdir -p "${OUT_DIR}"

####################################
# Collect grid files (simple & robust)
####################################
mapfile -t GRID_FILES < <(
    find "${GRID_DIR}" -maxdepth 1 -type f -name "*.h5ad" | sort
)

echo "Found ${#GRID_FILES[@]} grid files"

####################################
# Select file for this task
####################################
TASK_ID=${SLURM_ARRAY_TASK_ID:-0}

GRID_FILE="${GRID_FILES[$TASK_ID]}"

if [[ -z "${GRID_FILE}" ]]; then
    echo "No grid file for task ${TASK_ID}, exiting."
    exit 0
fi

GRID_BASENAME=$(basename "${GRID_FILE}" .h5ad)

echo "Running grid: ${GRID_BASENAME}"
echo "Input file:  ${GRID_FILE}"

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
# Run Python (logic lives here)
####################################
python << EOF
import os
import scanpy as sc
import stlearn as st

grid_file = "${GRID_FILE}"
out_dir = "${OUT_DIR}"
basename = os.path.basename(grid_file).replace(".h5ad", "")
out_path = os.path.join(out_dir, f"Sig_{basename}.h5ad")

print("Grid file:", grid_file)
print("Output:", out_path)

if os.path.exists(out_path):
    print("Output exists, skipping.")
    raise SystemExit(0)

grid = sc.read_h5ad(grid_file)
print("Grid shape:", grid.shape)

st.tl.cci.run_cci(
    grid,
    use_label="${USE_LABEL}",
    min_spots=${MIN_SPOTS},
    spot_mixtures=False,
    cell_prop_cutoff=${CELL_PROP_CUTOFF},
    sig_spots=True,
    n_perms=${N_PERMS},
    n_cpus=None
)

grid.write(out_path)
print("Saved:", out_path)
EOF

echo "Finished task ${TASK_ID}"
