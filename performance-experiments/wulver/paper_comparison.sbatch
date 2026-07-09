#!/bin/bash
#=============================================================
# Slurm SINGLE job: run the whole PAPER-COMPARISON suite (ac, ac-s2m2, cmg-v,
# cmg-k-elim; reps=3) sequentially in one allocation, then summarize into
# paper_comparison.{csv,md} + coverage.txt.
#
# Convenient for smoke/medium scale or a modest paper run on a big-memory node.
# For a full paper-scale sweep prefer the array path (faster, tiered memory):
#   CVK_REPS=3 ./run_paper_comparison.sh submit --scale paper --account <acct>
#
#   sbatch --export=ALL,CVK_SCALE=medium paper_comparison.sbatch
#   sbatch --mem=200G --time=72:00:00 --export=ALL,CVK_SCALE=paper paper_comparison.sbatch
#
# EDIT BEFORE FIRST USE: --account.
#=============================================================
#SBATCH --job-name=paper_comparison
#SBATCH --account=PI_UCID_HERE
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/paper_comparison_%j.out

set -u

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
    command -v julia >/dev/null 2>&1 || module load Julia/1.11.9 2>/dev/null || module load Julia 2>/dev/null || module load julia 2>/dev/null || true
fi
command -v julia >/dev/null 2>&1 || export PATH="$HOME/.juliaup/bin:$PATH"

# Keep the Julia depot off the small $HOME quota (see chol_vs_kcycle_array.sbatch);
# the Julia module clobbers JULIA_DEPOT_PATH, so set it AFTER module load.
export JULIA_DEPOT_PATH="${CVK_DEPOT:-/project/ikoutis/$USER/.julia}"

export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WULVER_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
EXP_DIR="$(dirname "$WULVER_DIR")"

"$EXP_DIR/run_paper_comparison.sh" run --scale "${CVK_SCALE:-medium}"
"$EXP_DIR/run_paper_comparison.sh" summarize
