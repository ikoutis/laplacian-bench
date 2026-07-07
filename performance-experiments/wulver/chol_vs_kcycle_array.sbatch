#!/bin/bash
#=============================================================
# Slurm ARRAY job for the chol-vs-kcycle benchmark on NJIT Wulver.
# One array task = one line of a manifest file (family [+ size] + options),
# as emitted by:  ../run_chol_vs_kcycle.sh --scale <S> --emit-manifest
#
# Submit from this directory (after ./setup_wulver.sh):
#   sbatch --array=1-$(wc -l < manifest.paper.small.txt) \
#          --export=ALL,CVK_MANIFEST=$PWD/manifest.paper.small.txt \
#          chol_vs_kcycle_array.sbatch
#   sbatch --array=1-$(wc -l < manifest.paper.large.txt) \
#          --mem=200G --time=72:00:00 \
#          --export=ALL,CVK_MANIFEST=$PWD/manifest.paper.large.txt \
#          chol_vs_kcycle_array.sbatch
#
# The header below is the SMALL tier (chimeras up to 1e6, medium grids).
# The LARGE tier (1e7 chimeras, nnz~2e8 grids, spe16m) is submitted with the
# --mem/--time overrides shown above; command-line flags override #SBATCH.
#
# EDIT BEFORE FIRST USE: --account (your PI's account UCID; list yours with
#   sacctmgr show associations user=$USER format=account%30 )
#=============================================================
#SBATCH --job-name=chol_vs_kcycle
#SBATCH --account=PI_UCID_HERE
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=logs/cvk_%A_%a.out

set -u

# Wulver module environment (module purge in a job requires re-loading wulver).
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
    # Pin Julia/1.11.9 (see setup_wulver.sh) so batch jobs match the setup env
    # and don't silently pick up the 1.12.6 default.
    command -v julia >/dev/null 2>&1 || module load Julia/1.11.9 2>/dev/null || module load Julia 2>/dev/null || module load julia 2>/dev/null || true
fi
command -v julia >/dev/null 2>&1 || export PATH="$HOME/.juliaup/bin:$PATH"

# Keep the Julia depot off the small $HOME quota. Set UNCONDITIONALLY and AFTER
# the module load: Wulver's EasyBuild Julia module (and --export=ALL forwarding)
# sets JULIA_DEPOT_PATH itself (to a value like ":"), so a ${JULIA_DEPOT_PATH:-...}
# fallback would keep that and Julia would miss our packages. Override via CVK_DEPOT.
export JULIA_DEPOT_PATH="${CVK_DEPOT:-/project/ikoutis/$USER/.julia}"

# All four solvers are single-threaded; pin the BLAS/OpenMP pools so timings
# are comparable and tasks stay inside their --cpus-per-task allocation.
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# When submitted from this directory (the documented flow), SLURM_SUBMIT_DIR
# is .../performance-experiments/wulver; fall back to the script location.
WULVER_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
EXP_DIR="$(dirname "$WULVER_DIR")"
ROOT="$(dirname "$EXP_DIR")"

MANIFEST="${CVK_MANIFEST:?set CVK_MANIFEST=/path/to/manifest.txt via --export}"
TASK_ARGS="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")"
if [[ -z "$TASK_ARGS" ]]; then
    echo "no task at line ${SLURM_ARRAY_TASK_ID} of $MANIFEST" >&2
    exit 1
fi

echo "host:   $(hostname)"
echo "depot:  ${JULIA_DEPOT_PATH:-<unset -> \$HOME/.julia>}"
echo "task:   ${SLURM_ARRAY_TASK_ID}: chol_vs_kcycle.jl $TASK_ARGS"
echo "julia:  $(command -v julia) ($(julia --version))"

cd "$EXP_DIR"
# shellcheck disable=SC2086  # TASK_ARGS is intentionally word-split
exec julia --project="$ROOT" chol_vs_kcycle.jl $TASK_ARGS
