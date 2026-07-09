#!/usr/bin/env bash
#
# Public paper-comparison run: ApproxChol (the paper's two solvers) vs CMG
# (legacy V-cycle + degree-1/2-elimination K-cycle) on the full benchmark,
# reproducing arXiv 2303.00709 (Gao-Kyng-Spielman) with CMG columns added.
#
# This is a THIN wrapper over the existing chol-vs-kcycle flow that pins the
# four public-table solver columns and the paper reps, so the run and the
# make_paper_tables.jl summary can't drift apart:
#
#   solvers = ac, ac-s2m2, cmg-v, cmg-k-elim      (--solvers, fixed here)
#   reps    = 3                                   (CVK_REPS, overridable)
#
# Usage (run on the Wulver LOGIN node after setup_wulver.sh has staged data):
#
#   ./run_paper_comparison.sh submit  [--scale paper] [--account ACCT] [extra submit args...]
#   ./run_paper_comparison.sh run     [--scale smoke]  [extra run args...]   # local/sequential
#   ./run_paper_comparison.sh summarize                                       # after jobs finish
#
#   submit    -> emit manifests + sbatch the Slurm array jobs (Wulver).
#   run       -> run every family sequentially on THIS machine (laptop smoke).
#   summarize -> aggregate the .jld2 results into paper_comparison.{csv,md} + coverage.txt.
#
# Env overrides: CVK_REPS (default 3), CVK_SEED (default 1), CVK_JULIA,
#                CVK_SCALE, CVK_DEPOT (Julia depot for setup), PAPER_SOLVERS.
#
# See wulver/README-paper-comparison.md for the full runbook (incl. the manual
# ipmMat.zip staging step and the intentionally-skipped SPE family).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$SCRIPT_DIR")"

# The four public-table columns (user decision). Override with PAPER_SOLVERS if
# you must, but the default IS the public table.
SOLVERS="${PAPER_SOLVERS:-ac,ac-s2m2,cmg-v,cmg-k-elim}"
REPS="${CVK_REPS:-3}"
SEED="${CVK_SEED:-1}"
JULIA_CMD="${CVK_JULIA:-julia}"

sub="${1:-}"
[[ $# -gt 0 ]] && shift || true

case "$sub" in
    submit)
        # Pass the pinned solvers + reps through the existing submit path.
        # CVK_EXTRA is threaded verbatim into every chol_vs_kcycle.jl task by
        # run_chol_vs_kcycle.sh; CVK_REPS overrides the per-scale default.
        echo "paper-comparison SUBMIT: solvers=$SOLVERS reps=$REPS seed=$SEED"
        CVK_EXTRA="--solvers $SOLVERS" CVK_REPS="$REPS" CVK_SEED="$SEED" \
            "$SCRIPT_DIR/wulver/submit_wulver.sh" "$@"
        ;;
    run)
        # Sequential local run (laptop/smoke). Defaults to smoke scale.
        echo "paper-comparison RUN (local): solvers=$SOLVERS reps=$REPS seed=$SEED"
        CVK_EXTRA="--solvers $SOLVERS" CVK_REPS="$REPS" CVK_SEED="$SEED" \
            "$SCRIPT_DIR/run_chol_vs_kcycle.sh" --run "$@"
        ;;
    summarize)
        # Summarize ONLY this run's result files, matched by the seed/reps
        # signature, so stale files from other runs sharing the results
        # directory (older smoke/reps1 dev runs, committed sample results) do not
        # pollute the coverage report. Pass explicit files/dirs to override.
        RESULTS_DIR="$ROOT/performance-analyses/chol-vs-kcycle"
        echo "paper-comparison SUMMARIZE: solvers=$SOLVERS (seed=$SEED reps=$REPS)"
        if [[ $# -gt 0 ]]; then
            (cd "$SCRIPT_DIR" && "$JULIA_CMD" --project="$ROOT" make_paper_tables.jl \
                --solvers "$SOLVERS" "$@")
        else
            shopt -s nullglob
            FILES=( "$RESULTS_DIR"/*.seed${SEED}.reps${REPS}.jld2 )
            shopt -u nullglob
            if [[ ${#FILES[@]} -eq 0 ]]; then
                echo "no result files matching *.seed${SEED}.reps${REPS}.jld2 in $RESULTS_DIR" >&2
                echo "(check --reps/--seed, or pass files/dirs explicitly)" >&2
                exit 1
            fi
            echo "summarizing ${#FILES[@]} file(s) matching *.seed${SEED}.reps${REPS}.jld2"
            (cd "$SCRIPT_DIR" && "$JULIA_CMD" --project="$ROOT" make_paper_tables.jl \
                --solvers "$SOLVERS" "${FILES[@]}")
        fi
        ;;
    selftest)
        # Fast, data-free gate for the summarizer (laptop CI).
        (cd "$SCRIPT_DIR" && "$JULIA_CMD" --project="$ROOT" make_paper_tables.jl --selftest)
        ;;
    -h|--help|"")
        sed -n '2,40p' "$0"
        ;;
    *)
        echo "unknown subcommand: $sub (want submit|run|summarize|selftest)" >&2
        exit 1
        ;;
esac
