#!/usr/bin/env bash
#
# Load Julia + the project's package depot for an INTERACTIVE Wulver session,
# so you don't retype the module/depot setup each time. SOURCE it (don't run it):
#
#     source performance-experiments/wulver/env_wulver.sh
#     julia --project=. -e '...'
#
# Works on a login node or (after `srun ... --pty bash`) a compute node. It only
# sets up the environment — it does NOT instantiate/precompile packages; run
# setup_wulver.sh once on the login node first for that. Mirrors what
# setup_wulver.sh (module + depot) and the sbatch scripts (thread pins) do.
#
# Overrides (export before sourcing):
#   CVK_DEPOT           Julia depot location (default /project/ikoutis/$USER/.julia)
#   JULIA_PKG_OFFLINE   set to false on a login node if you need to install pkgs

# 1. Wulver module environment; load Julia only if it isn't already on PATH.
if command -v module >/dev/null 2>&1; then
    module purge || true
    module load wulver || true
    command -v julia >/dev/null 2>&1 || \
        module load Julia/1.11.9 2>/dev/null || \
        module load Julia 2>/dev/null || \
        module load julia 2>/dev/null || true
fi
command -v julia >/dev/null 2>&1 || [[ ! -x "$HOME/.juliaup/bin/julia" ]] || \
    export PATH="$HOME/.juliaup/bin:$PATH"

# 2. Package depot (must match what setup_wulver.sh used, else Julia won't find
#    the instantiated packages). The EasyBuild Julia module clobbers
#    JULIA_DEPOT_PATH, so set it AFTER the module load above.
export JULIA_DEPOT_PATH="${CVK_DEPOT:-/project/ikoutis/$USER/.julia}"

# 3. Compute nodes have no internet; keep Julia from reaching the registry.
#    Override with JULIA_PKG_OFFLINE=false on the login node for installs.
export JULIA_PKG_OFFLINE="${JULIA_PKG_OFFLINE:-true}"

# 4. Single-threaded, matching the benchmark sbatch scripts (fair timings).
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "env_wulver: julia=$(command -v julia || echo MISSING) depot=$JULIA_DEPOT_PATH offline=$JULIA_PKG_OFFLINE"
