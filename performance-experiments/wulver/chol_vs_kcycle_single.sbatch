#!/bin/bash
#=============================================================
# Slurm SINGLE job: run the whole chol-vs-kcycle suite sequentially in one
# allocation. Convenient for smoke/medium scale; use the array job
# (chol_vs_kcycle_array.sbatch) for paper scale.
#
#   sbatch --export=ALL,CVK_SCALE=medium chol_vs_kcycle_single.sbatch
#
# EDIT BEFORE FIRST USE: --account (as in chol_vs_kcycle_array.sbatch).
#=============================================================
#SBATCH --job-name=chol_vs_kcycle_seq
#SBATCH --account=PI_UCID_HERE
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/cvk_seq_%j.out

set -u

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
    # Pin Julia/1.11.9 (see setup_wulver.sh) so batch jobs match the setup env
    # and don't silently pick up the 1.12.6 default.
    command -v julia >/dev/null 2>&1 || module load Julia/1.11.9 2>/dev/null || module load Julia 2>/dev/null || module load julia 2>/dev/null || true
fi
command -v julia >/dev/null 2>&1 || export PATH="$HOME/.juliaup/bin:$PATH"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WULVER_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
EXP_DIR="$(dirname "$WULVER_DIR")"

exec "$EXP_DIR/run_chol_vs_kcycle.sh" --scale "${CVK_SCALE:-medium}" --run
