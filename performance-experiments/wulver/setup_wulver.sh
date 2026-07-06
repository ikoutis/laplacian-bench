#!/usr/bin/env bash
#
# One-time setup on NJIT Wulver (run on a LOGIN node — compute nodes have no
# internet access, and both the Julia packages and the benchmark matrices are
# fetched here).
#
#   cd performance-experiments/wulver && ./setup_wulver.sh [--scale paper]
#
# Steps:
#   1. load the Wulver module environment,
#   2. make Julia available (site module if present, else juliaup in $HOME,
#      pinned to the 1.10 LTS channel to match the pinned environment),
#   3. instantiate the Julia project (adds CombinatorialMultigrid from GitHub
#      via ../..//setup.jl on the first run, or plain Pkg.instantiate() once a
#      Manifest.toml is committed),
#   4. prefetch the downloadable benchmark matrices,
#   5. create the logs/ directory the sbatch scripts write into.
#
# If your $HOME quota is tight, point the Julia depot somewhere roomier
# BEFORE running this, e.g.:  export JULIA_DEPOT_PATH=/project/<PI>/$USER/.julia

set -euo pipefail

SCALE="paper"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --scale) SCALE="$2"; shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP_DIR="$(dirname "$SCRIPT_DIR")"
ROOT="$(dirname "$EXP_DIR")"

# 1. Wulver module environment (harmless elsewhere).
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
fi

# 2. Julia: site module, existing binary, or juliaup.
if ! command -v julia >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    # Pin the Julia module: Project.toml requires >=1.10 and Wulver has no 1.10
    # module (only 1.9.3, 1.11.x, 1.12.x). 1.11.9 is the newest stable release
    # the deps precompile cleanly on; the bare default is 1.12.6, which is
    # riskier for older packages like Laplacians. Fall back if 1.11.9 is absent.
    module load Julia/1.11.9 2>/dev/null || module load Julia 2>/dev/null || module load julia 2>/dev/null || true
fi

# Keep the Julia depot off the small $HOME quota. Set UNCONDITIONALLY and AFTER
# the module load: Wulver's EasyBuild Julia module sets JULIA_DEPOT_PATH itself
# (to a value like ":"), so a ${JULIA_DEPOT_PATH:-...} fallback would keep the
# module's value and miss our packages. Override the location via CVK_DEPOT.
export JULIA_DEPOT_PATH="${CVK_DEPOT:-/project/ikoutis/$USER/.julia}"
echo "using JULIA_DEPOT_PATH=$JULIA_DEPOT_PATH"

if ! command -v julia >/dev/null 2>&1 && [[ -x "$HOME/.juliaup/bin/julia" ]]; then
    export PATH="$HOME/.juliaup/bin:$PATH"
fi
if ! command -v julia >/dev/null 2>&1; then
    echo "no julia found — installing juliaup (channel 1.10) into \$HOME"
    curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel 1.10
    export PATH="$HOME/.juliaup/bin:$PATH"
fi
julia --version

# 3. Julia project.
cd "$ROOT"
if [[ -f Manifest.toml ]]; then
    echo "Manifest.toml present — instantiating the pinned environment"
    julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
else
    echo "no Manifest.toml — running setup.jl (adds CombinatorialMultigrid from GitHub)"
    julia --project=. setup.jl
    echo "consider committing the generated Manifest.toml to pin this environment"
fi

# 4. Benchmark data (SuiteSparse + IPM downloads; reports missing SPE files).
cd "$EXP_DIR"
julia --project="$ROOT" download_data.jl --scale "$SCALE" || {
    echo "WARNING: some downloads failed; the affected instances will be skipped" >&2
}

# 5. Log directory for sbatch output.
mkdir -p "$SCRIPT_DIR/logs"

echo
echo "setup complete. next:"
echo "  1. edit chol_vs_kcycle_array.sbatch: set --account to your PI account"
echo "     (find it with: sacctmgr show associations user=\$USER format=account%30)"
echo "  2. ./submit_wulver.sh --scale $SCALE"
