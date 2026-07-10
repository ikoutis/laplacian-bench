#!/usr/bin/env bash
#
# Emit the task manifests for the requested scale and submit the Slurm array
# job(s) on Wulver — small tier with the sbatch-file defaults, large tier with
# the big-memory overrides.
#
#   ./submit_wulver.sh [--scale paper] [--reps R] [--seed S]
#                      [--account ACCT] [--sbatch-extra "..."] [--dry-run]
#
# --account overrides the placeholder in the .sbatch file for this submission
# only (sbatch command-line flags beat #SBATCH headers).

set -euo pipefail

SCALE="paper"
REPS=""
SEED=""
ACCOUNT=""
SBATCH_EXTRA=""
DRY_RUN=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --scale) SCALE="$2"; shift 2 ;;
        --reps) REPS="$2"; shift 2 ;;
        --seed) SEED="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --sbatch-extra) SBATCH_EXTRA="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) sed -n '2,13p' "$0"; exit 0 ;;
        *) echo "unknown argument: $1" >&2; exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXP_DIR="$(dirname "$SCRIPT_DIR")"

EMIT_ARGS=(--scale "$SCALE" --emit-manifest)
[[ -n "$REPS" ]] && EMIT_ARGS+=(--reps "$REPS")
[[ -n "$SEED" ]] && EMIT_ARGS+=(--seed "$SEED")
"$EXP_DIR/run_chol_vs_kcycle.sh" "${EMIT_ARGS[@]}"

mkdir -p "$SCRIPT_DIR/logs"

SBATCH_ARGS=()
[[ -n "$ACCOUNT" ]] && SBATCH_ARGS+=(--account "$ACCOUNT")
# shellcheck disable=SC2206  # intentional word-splitting of extra flags
[[ -n "$SBATCH_EXTRA" ]] && SBATCH_ARGS+=($SBATCH_EXTRA)

submit() {
    echo "+ sbatch $*"
    if [[ "$DRY_RUN" == "0" ]]; then
        (cd "$SCRIPT_DIR" && sbatch "$@")
    fi
}

SMALL_MF="$SCRIPT_DIR/manifest.$SCALE.small.txt"
LARGE_MF="$SCRIPT_DIR/manifest.$SCALE.large.txt"

if [[ -s "$SMALL_MF" ]]; then
    N=$(wc -l < "$SMALL_MF")
    submit --array="1-$N" --export=ALL,CVK_MANIFEST="$SMALL_MF" \
        "${SBATCH_ARGS[@]+"${SBATCH_ARGS[@]}"}" chol_vs_kcycle_array.sbatch
fi
if [[ -s "$LARGE_MF" ]]; then
    N=$(wc -l < "$LARGE_MF")
    # 48h fits the chunked 1e6/1e7 chimeras (each element runs a few whole
    # instances). A from-scratch run whose large tier also holds the biggest
    # un-chunked fixed graphs (spielmanIPM k500/k600, nnz~2e8 grids) can add
    # more time via --sbatch-extra "--time=72:00:00".
    submit --array="1-$N" --mem=200G --time=48:00:00 --export=ALL,CVK_MANIFEST="$LARGE_MF" \
        "${SBATCH_ARGS[@]+"${SBATCH_ARGS[@]}"}" chol_vs_kcycle_array.sbatch
fi

echo "monitor with: squeue -u \$USER ; results land in performance-analyses/chol-vs-kcycle/"
