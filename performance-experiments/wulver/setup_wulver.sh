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

# 2b. Fetch CombinatorialMultigrid.jl into the repo (pure git — safe here on the
#     login node). setup.jl then develops this local checkout instead of doing a
#     run-time GitHub fetch, so compute nodes need no internet. Override the
#     branch/tag with CMG_FETCH_REV (default main). Guarded so a transient git
#     failure doesn't abort setup when a usable checkout already exists.
"$SCRIPT_DIR/fetch_cmg.sh" --rev "${CMG_FETCH_REV:-main}" \
    || echo "WARNING: fetch_cmg.sh failed; continuing with any existing checkout" >&2

# 3. Julia project. Precompilation is deferred (SETUP_SKIP_PRECOMPILE=1): it
#    happens on first `using` on a compute node, so this step stays light — some
#    clusters (e.g. Wulver) restrict heavy compute on the login node. If even the
#    instantiate below fails here, rerun it on a compute node (see the note at the
#    end); the in-repo CMG checkout makes it fully offline.
cd "$ROOT"
export SETUP_SKIP_PRECOMPILE=1
SETUP_JULIA_FAILED=0
if [[ -f Manifest.toml ]]; then
    echo "Manifest.toml present — instantiating (precompile deferred to compute)"
    julia --project=. -e 'using Pkg; Pkg.instantiate()' || SETUP_JULIA_FAILED=1
else
    echo "no Manifest.toml — running setup.jl (develops the in-repo CMG checkout)"
    julia --project=. setup.jl || SETUP_JULIA_FAILED=1
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
if [[ "$SETUP_JULIA_FAILED" == "1" ]]; then
    cat <<NOTE
NOTE: the Julia instantiate step did not finish on this (login) node.
 * If it downloaded everything and only the compile/precompile choked, rerun it
   on a COMPUTE node — the in-repo CombinatorialMultigrid.jl checkout makes that
   fully offline:
       srun --account=<acct> --partition=general --qos=standard \\
            --cpus-per-task=4 --mem=32G --time=01:00:00 --pty bash
       cd "$ROOT"
       export JULIA_DEPOT_PATH=<your depot, e.g. /project/<PI>/\$USER/.julia>
       JULIA_PKG_OFFLINE=true julia --project=. setup.jl
 * If it failed while DOWNLOADING (couldn't reach the network), fix connectivity
   and rerun THIS script on the login node first — an offline compute node can
   never fetch CombinatorialMultigrid's deps (BenchmarkTools, LDLFactorizations,
   artifacts) if they aren't already in your depot.

NOTE
fi
echo "setup complete. next:"
echo "  1. edit chol_vs_kcycle_array.sbatch: set --account to your PI account"
echo "     (find it with: sacctmgr show associations user=\$USER format=account%30)"
echo "  2. ./submit_wulver.sh --scale $SCALE"
echo
echo "  or, for the Spielman degree-1/2 elimination comparison alone"
echo "  (ac, ac-s2m2, cmg-k, cmg-k-elim with per-matrix core%):"
echo "    edit spielman_elim.sbatch (--account), then:"
echo "    sbatch --export=ALL,SPL_SCALE=$SCALE spielman_elim.sbatch"
