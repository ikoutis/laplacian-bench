#!/usr/bin/env bash
#
# Superscript for the approx-Cholesky vs CMG K-cycle benchmark.
#
# Runs every benchmark family (chimera families once per size) through
# performance-experiments/chol_vs_kcycle.jl, either sequentially on this
# machine (--run, the default) or by emitting Slurm-array manifests
# (--emit-manifest) consumed by wulver/chol_vs_kcycle_array.sbatch.
#
# Usage:
#   ./run_chol_vs_kcycle.sh [--scale smoke|medium|paper] [--reps R] [--seed S]
#                           [--run | --emit-manifest] [--julia CMD]
#                           [--extra "ARGS..."] [--dry-run]
#
# Defaults: --scale smoke; reps by scale (paper 3, medium 5, smoke 2); seed 1.
# Environment overrides: CVK_SCALE, CVK_REPS, CVK_SEED, CVK_JULIA, CVK_EXTRA.
#
# Examples:
#   ./run_chol_vs_kcycle.sh                                   # quick smoke pass
#   ./run_chol_vs_kcycle.sh --scale medium --reps 5 --run
#   ./run_chol_vs_kcycle.sh --scale paper --emit-manifest     # for Wulver
#   CVK_EXTRA="--solvers ac,cmg-k" ./run_chol_vs_kcycle.sh    # subset of columns

set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

SCALE="${CVK_SCALE:-smoke}"
REPS="${CVK_REPS:-}"
SEED="${CVK_SEED:-1}"
JULIA_CMD="${CVK_JULIA:-julia}"
EXTRA="${CVK_EXTRA:-}"
MODE="run"
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scale)  SCALE="$2"; shift 2 ;;
        --reps)   REPS="$2"; shift 2 ;;
        --seed)   SEED="$2"; shift 2 ;;
        --julia)  JULIA_CMD="$2"; shift 2 ;;
        --extra)  EXTRA="$2"; shift 2 ;;
        --run)    MODE="run"; shift ;;
        --emit-manifest) MODE="manifest"; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,25p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

case "$SCALE" in
    smoke)  DEFAULT_REPS=2; CHIMERA_SIZES=(1e4) ;;
    medium) DEFAULT_REPS=5; CHIMERA_SIZES=(1e4 1e5) ;;
    paper)  DEFAULT_REPS=3; CHIMERA_SIZES=(1e4 1e5 1e6 1e7) ;;
    *) echo "bad --scale $SCALE (want smoke|medium|paper)" >&2; exit 1 ;;
esac
REPS="${REPS:-$DEFAULT_REPS}"

COMMON="--scale $SCALE --reps $REPS --seed $SEED"
[[ -n "$EXTRA" ]] && COMMON="$COMMON $EXTRA"

# Families with a fixed instance sweep (one task each).
FIXED_FAMILIES=(uniform_grid aniso wgrid checkered sachdeva_star suitesparse spe chimeraIPM spielmanIPM)
# Chimera families run once per size (one task per family x size).
CHIMERA_FAMILIES=(uni_chimera uni_bndry_chimera wted_chimera wted_bndry_chimera)
# At paper scale these need the big-memory tier (nnz~2e8 grids, spe16m, 1e7 chimeras).
BIG_FIXED="aniso wgrid checkered uniform_grid spe"

SMALL_TASKS=()
LARGE_TASKS=()

for fam in "${FIXED_FAMILIES[@]}"; do
    task="$fam $COMMON"
    if [[ "$SCALE" == "paper" && " $BIG_FIXED " == *" $fam "* ]]; then
        LARGE_TASKS+=("$task")
    else
        SMALL_TASKS+=("$task")
    fi
done
for fam in "${CHIMERA_FAMILIES[@]}"; do
    for n in "${CHIMERA_SIZES[@]}"; do
        task="$fam --n $n $COMMON"
        if [[ "$SCALE" == "paper" && "$n" == "1e7" ]]; then
            LARGE_TASKS+=("$task")
        else
            SMALL_TASKS+=("$task")
        fi
    done
done

echo "scale=$SCALE reps=$REPS seed=$SEED mode=$MODE"
echo "tasks: ${#SMALL_TASKS[@]} small + ${#LARGE_TASKS[@]} large"

if [[ "$MODE" == "manifest" ]]; then
    OUTDIR="$SCRIPT_DIR/wulver"
    mkdir -p "$OUTDIR"
    SMALL_MF="$OUTDIR/manifest.$SCALE.small.txt"
    LARGE_MF="$OUTDIR/manifest.$SCALE.large.txt"
    : > "$SMALL_MF"
    : > "$LARGE_MF"
    for t in "${SMALL_TASKS[@]+"${SMALL_TASKS[@]}"}"; do echo "$t" >> "$SMALL_MF"; done
    for t in "${LARGE_TASKS[@]+"${LARGE_TASKS[@]}"}"; do echo "$t" >> "$LARGE_MF"; done
    echo "wrote $SMALL_MF ($(wc -l < "$SMALL_MF") tasks)"
    echo "wrote $LARGE_MF ($(wc -l < "$LARGE_MF") tasks)"
    echo
    echo "submit on Wulver (edit the account in the sbatch file first):"
    echo "  cd $SCRIPT_DIR/wulver"
    if [[ -s "$SMALL_MF" ]]; then
        echo "  sbatch --array=1-$(wc -l < "$SMALL_MF") --export=ALL,CVK_MANIFEST=$SMALL_MF chol_vs_kcycle_array.sbatch"
    fi
    if [[ -s "$LARGE_MF" ]]; then
        echo "  sbatch --array=1-$(wc -l < "$LARGE_MF") --mem=200G --time=72:00:00 --export=ALL,CVK_MANIFEST=$LARGE_MF chol_vs_kcycle_array.sbatch"
    fi
    exit 0
fi

# --run mode: sequential execution, small tier first.
FAILED=0
run_task() {
    local task="$1"
    echo
    echo "=================================================================="
    echo ">>> chol_vs_kcycle.jl $task"
    echo "=================================================================="
    if [[ "$DRY_RUN" == "1" ]]; then
        return 0
    fi
    # shellcheck disable=SC2086  # task is intentionally word-split into args
    (cd "$SCRIPT_DIR" && $JULIA_CMD --project="$ROOT" chol_vs_kcycle.jl $task)
    local rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "!!! task failed (rc=$rc): $task" >&2
        FAILED=$((FAILED + 1))
    fi
    return 0
}

for t in "${SMALL_TASKS[@]+"${SMALL_TASKS[@]}"}"; do run_task "$t"; done
for t in "${LARGE_TASKS[@]+"${LARGE_TASKS[@]}"}"; do run_task "$t"; done

echo
if [[ $FAILED -gt 0 ]]; then
    echo "$FAILED task(s) failed" >&2
    exit 1
fi
echo "all tasks completed"
