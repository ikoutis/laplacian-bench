# Wulver runbook v2: measure the fast elimination on the Spielman graphs

## Context

CMG.jl `main` (`f0f5e40`) has the array-based `eliminate_deg12` (PR #3;
targets the 10.8 s k300 elimination build) and honest-failure breakdown
semantics (the guarded-refinement salvage was reverted in PR #5).
laplacian-bench `main` (`80c4caf` and later — seeing this v2 text after a pull
confirms you're current) has the anonymous-environment Project.toml fix, so
setup steps exit 0. Goal: re-run the Spielman comparison on a compute node and
compare `cmg-k-elim`'s `build_s` against the previous run
(0.300 / 2.490 / 10.772 s at k100/k200/k300) and against `ac`
(0.208 / 0.845 / 3.031 s).

Landmines, unchanged: login = git only (Julia OOMs on precompile), compute =
offline, depot exported **after** `module load` (the Julia module clobbers
`JULIA_DEPOT_PATH`). New simplification: the bench env `develop`s the in-repo
CMG checkout, so after a `git pull` of that checkout **no Pkg operation is
needed** — the first `using` on the compute node recompiles the changed source.

Run the phases in order; don't proceed past a failed check.

## Phase A — login node (git only, no Julia)

```bash
cd /project/ikoutis/github/laplacian-bench
git status --short                                        # if dirty (result files):
git stash push -m "local results before fast-elim run"    #   only if dirty
git checkout main
git pull --ff-only origin main

performance-experiments/wulver/fetch_cmg.sh               # in-repo CMG -> main
git -C CombinatorialMultigrid.jl log --oneline -1         # EXPECT: f0f5e40 Merge pull request #5 ...
```

## Phase B — grab a compute node

```bash
srun --account=ikoutis --partition=general --qos=standard \
     --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G --time=04:00:00 --pty bash
```

(64G covers the medium sweep k100–k300; k500/k600 are 126M/217M nnz — Phase D.)

## Phase C — on the compute node

```bash
cd /project/ikoutis/github/laplacian-bench
module purge; module load wulver
command -v julia >/dev/null || module load Julia/1.11.9
export JULIA_DEPOT_PATH=/project/ikoutis/ikoutis/.julia   # AFTER module load
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "$JULIA_DEPOT_PATH"                                  # must not be ":" or empty

# gate: first `using` recompiles CMG from the updated in-repo checkout;
# compact_adjacency! exists only in the new fast elimination
julia --project=. -e '
  using CombinatorialMultigrid
  @show pathof(CombinatorialMultigrid)
  @show isdefined(CombinatorialMultigrid, :compact_adjacency!)'
# pathof must be under .../laplacian-bench/CombinatorialMultigrid.jl/
# isdefined must be true — do not benchmark before both hold

# measurement run (Spielman_k100/k200/k300)
julia --project=. performance-experiments/spielman_elim_demo.jl --scale medium --offline \
  | tee performance-experiments/wulver/spielman_medium_fastelim.out
```

Failure handling:
- `pathof` points into `$JULIA_DEPOT_PATH/packages/...` → the Manifest lost the
  develop entry; run `JULIA_PKG_OFFLINE=true julia --project=. setup.jl`
  (auto-develops the in-repo checkout, exits 0) and re-check.
- `isdefined(... :compact_adjacency!)` is false → the in-repo checkout wasn't
  updated; redo Phase A's `fetch_cmg.sh` on login.
- "Package ... not installed" → depot wrong; redo the module-load-then-export
  order above.
- Precompile hang / `✗ Pkg` → `rm -rf $JULIA_DEPOT_PATH/compiled/v1.11/Pkg`,
  retry.

## Phase D — optional paper scale (k100..k600), after C succeeds

```bash
cd performance-experiments/wulver && mkdir -p logs
sbatch --export=ALL,SPL_SCALE=paper,SPL_MAXITS=2000 --mem=200G --time=48:00:00 spielman_elim.sbatch
```

## What to look for

- `cmg-k-elim` `build_s` vs the previous 0.300 / 2.490 / 10.772 s — target ≥5×
  lower (toward or below ac's 0.208 / 0.845 / 3.031 s). `solve_s` should stay
  ~0.007 / 0.055 / 0.19 s and `relres` at 1e-10…1e-14.
- `core` must still be exactly 100 / 200 / 300 (algorithm unchanged, only its
  bookkeeping).
- Expected and fine: k300 `cmg-k-elim` prints `converged=false` at relres
  ~1.2e-8 (honest breakdown semantics, by choice); the table's `relres` column
  is the ground truth.
- Comparison caveat: a different compute node than the last run (n0056) can
  shift timings by tens of percent; the build-time *ratio* to `ac` within the
  same run is the robust number.
- Decision afterwards: if `cmg-k-elim` build ≈ ac build, the remaining
  optimization ideas (leaf sweep first, pooled records, chain-walking) are
  optional polish; otherwise pick them up next.
