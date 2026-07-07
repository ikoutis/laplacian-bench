#!/bin/bash
#=============================================================
# Slurm job: run spielman_elim_demo.jl on a Wulver compute node — the
# degree-1/2 elimination comparison (ac, ac-s2m2, cmg-k, cmg-k-elim) on the
# Spielman IPM Laplacians, with the per-matrix core% report.
#
# Prereq: run wulver/setup_wulver.sh on a LOGIN node first — it instantiates the
# Julia project (via ../../setup.jl) and prefetches the sk<k>i<i> matrices into
# matrix-files/. Compute nodes have no internet, so this job runs --offline and
# uses only that cache. Make sure setup pointed at the CMG `eliminate` branch,
# e.g.  CMG_REV=main ./setup_wulver.sh --scale paper.
#
#   sbatch --export=ALL,SPL_SCALE=medium spielman_elim.sbatch
#   sbatch --export=ALL,SPL_SCALE=paper,SPL_MAXITS=2000 --mem=96G --time=24:00:00 spielman_elim.sbatch
#   sbatch --export=ALL,SPL_NAMES=sk100i1,sk300i5,sk600i11 spielman_elim.sbatch
#
# Knobs (all optional, via --export):
#   SPL_SCALE    smoke|medium|paper                       [medium]
#   SPL_NAMES    explicit sk<k>i<i> list (overrides SPL_SCALE)
#   SPL_LIMIT    cap number of matrices
#   SPL_SOLVERS  column subset  [ac,ac-s2m2,cmg-k,cmg-k-elim]
#   SPL_MAXITS   max outer iterations                     [1000]
#   SPL_TOL      residual tolerance                       [1e-8]
#   SPL_SEED     base seed                                [1]
#   SPL_EXTRA    extra raw args passed through verbatim
#
# EDIT BEFORE FIRST USE: set --account to your PI account
#   (find it with: sacctmgr show associations user=$USER format=account%30)
#=============================================================
#SBATCH --job-name=spielman_elim
#SBATCH --account=PI_UCID_HERE
#SBATCH --partition=general
#SBATCH --qos=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=logs/spielman_elim_%j.out

set -u

if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
    # Pin the same Julia the setup step used (see setup_wulver.sh) so the batch
    # job matches the instantiated environment and doesn't pick up the default.
    command -v julia >/dev/null 2>&1 || module load Julia/1.11.9 2>/dev/null || module load Julia 2>/dev/null || module load julia 2>/dev/null || true
fi
command -v julia >/dev/null 2>&1 || export PATH="$HOME/.juliaup/bin:$PATH"

# Single-threaded, matching the other Wulver jobs, so solver timings are
# comparable across columns and across runs.
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WULVER_DIR="${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}"
EXP_DIR="$(dirname "$WULVER_DIR")"
ROOT="$(dirname "$EXP_DIR")"

# Build the argument list from the SPL_* knobs.
ARGS=()
if [[ -n "${SPL_NAMES:-}" ]]; then
    ARGS+=(--names "$SPL_NAMES")
else
    ARGS+=(--scale "${SPL_SCALE:-medium}")
fi
[[ -n "${SPL_LIMIT:-}" ]]   && ARGS+=(--limit "$SPL_LIMIT")
[[ -n "${SPL_SOLVERS:-}" ]] && ARGS+=(--solvers "$SPL_SOLVERS")
[[ -n "${SPL_MAXITS:-}" ]]  && ARGS+=(--maxits "$SPL_MAXITS")
[[ -n "${SPL_TOL:-}" ]]     && ARGS+=(--tol "$SPL_TOL")
[[ -n "${SPL_SEED:-}" ]]    && ARGS+=(--seed "$SPL_SEED")
ARGS+=(--offline)
# shellcheck disable=SC2206  # SPL_EXTRA is intentionally word-split into args
[[ -n "${SPL_EXTRA:-}" ]]   && ARGS+=(${SPL_EXTRA})

echo "julia: $(command -v julia)  ($(julia --version 2>/dev/null))"
echo "args: ${ARGS[*]}"

exec julia --project="$ROOT" "$EXP_DIR/spielman_elim_demo.jl" "${ARGS[@]}"
