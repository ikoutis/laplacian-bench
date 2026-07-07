# Wulver runbook: test CMG degree-1/2 elimination on the Spielman graphs

## Context

Everything is merged: laplacian-bench `main` (`1b61792` and later — seeing this
file after a pull confirms you're current) has the elimination columns
(`cmg-k-elim`/`cmg-v-elim`), `spielman_elim_demo.jl`, the FlowIPM22 Spielman
loaders (`Spielman_k100..k600`), the Wulver depot fixes, and in-repo CMG
support; CombinatorialMultigrid `main` (`98fe870`) has the `eliminate` branch. The goal now is purely operational: on Wulver, get onto a compute node
and run the comparison (`ac`, `ac-s2m2`, `cmg-k`, `cmg-k-elim` + per-matrix
`core%`) on the real Spielman graphs, without tripping the known landmines:

- **login node**: has internet, but OOMs on Julia precompile → git + downloads only.
- **compute node**: no internet → all Julia (setup, precompile, run) happens here, offline.
- **depot**: Wulver's Julia module clobbers `JULIA_DEPOT_PATH` (sets it to `:`), so it
  must be exported **after** `module load`, unconditionally, to
  `/project/ikoutis/ikoutis/.julia` (the depot holding all previously installed deps).
- Paths: repo = `/project/ikoutis/github/laplacian-bench`.

This is an exact command sequence with checks between steps — run the phases in
order and don't proceed past a failed check.

## Phase A — login node (login02): git only, no Julia

```bash
cd /project/ikoutis/github/laplacian-bench

# A1. clean tree state (local .jld2 results block checkout otherwise)
git status --short                 # see what's dirty
git stash push -m "local results before elim run"   # only if dirty
git checkout main
git pull --ff-only origin main
ls performance-experiments/wulver/RUN-ELIM-ON-SPIELMAN.md   # EXPECT: exists (you're current)

# A2. update the in-repo CMG clone to main (pure git; tolerant of hiccups)
performance-experiments/wulver/fetch_cmg.sh
ls CombinatorialMultigrid.jl/src/elimination.jl    # MUST exist — abort if not

# A3. confirm the Spielman matrices are cached (prefetched earlier)
ls -lh matrix-files/Spielman_k*.mat matrix-files/Spielman_k*.mm 2>/dev/null
```

**Check A3:** for `--scale medium` you need `Spielman_k100/k200/k300`. If any are
missing, fetch them on login with plain HTTP (no Julia):
`wget -O matrix-files/Spielman_k<k>.mat http://sparse-files.engr.tamu.edu/mat/FlowIPM22/Spielman_k<k>.mat`

## Phase B — grab an interactive compute node

```bash
srun --account=ikoutis --partition=general --qos=standard \
     --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G --time=04:00:00 --pty bash
```

(64G covers k300 comfortably; k500/k600 are 126M/217M nnz and are NOT in the
medium sweep — they need a big-memory job, Phase D.)

## Phase C — on the compute node: environment, verify, run

```bash
cd /project/ikoutis/github/laplacian-bench

# C1. Julia + depot — ORDER MATTERS: depot export AFTER module load
module purge; module load wulver
command -v julia >/dev/null || module load Julia/1.11.9
export JULIA_DEPOT_PATH=/project/ikoutis/ikoutis/.julia
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
echo "$JULIA_DEPOT_PATH"; julia --version     # depot must NOT be ":" or empty

# C2. point the env at the in-repo CMG and precompile — fully offline
#     (setup.jl auto-develops ./CombinatorialMultigrid.jl; all deps are already
#      in the depot from earlier runs, so no network is needed or attempted)
JULIA_PKG_OFFLINE=true julia --project=. setup.jl

# C3. verify the NEW CMG is what's loaded — do not benchmark before this passes
julia --project=. -e '
  using CombinatorialMultigrid
  @show pathof(CombinatorialMultigrid)                       # must be under /project/ikoutis/github/laplacian-bench/CombinatorialMultigrid.jl/
  @show isdefined(CombinatorialMultigrid, :EliminatedHierarchy)   # must be true'

# C4. run the comparison (medium = Spielman_k100,k200,k300), save output
julia --project=. performance-experiments/spielman_elim_demo.jl --scale medium --offline \
  | tee performance-experiments/wulver/spielman_medium.out
```

Failure handling:
- C2 errors mentioning a package "not installed" → depot is wrong; recheck C1
  (the module clobbered it, or the export ran before `module load`).
- C2 hangs on "waiting for lock" or reprints `✗ Pkg` → stale cache from the
  earlier killed login precompile: `rm -rf $JULIA_DEPOT_PATH/compiled/v1.11/Pkg`
  and rerun C2.
- C3 `pathof` pointing into `$JULIA_DEPOT_PATH/packages/...` → the env still has
  the URL-pinned CMG; rerun C2 and confirm it logs "Using in-repo
  CombinatorialMultigrid checkout".
- Demo prints "No Spielman matrices available" → cache check A3 failed for the
  bare names; fetch on login as in A3, rerun (it is `--offline`, so it never
  downloads on compute).

## Phase D — optional: full paper sweep via Slurm batch (k100..k600)

Only after Phase C succeeds. `Spielman_k500/k600` are the large ones:

```bash
# from login or compute shell, in performance-experiments/wulver/
mkdir -p logs
sbatch --export=ALL,SPL_SCALE=paper,SPL_MAXITS=2000 --mem=200G --time=48:00:00 spielman_elim.sbatch
# or skip the giants first:  --export=ALL,SPL_SCALE=paper,SPL_LIMIT=4
tail -f logs/spielman_elim_<jobid>.out
```

The sbatch already handles module/depot/threads itself (account=ikoutis,
depot default `/project/ikoutis/$USER/.julia`, `--offline` hard-coded).

## Verification / what success looks like

- C3 prints the in-repo path and `true`.
- `spielman_medium.out` shows, per matrix, a `core= ... (x.xx% of n)` header —
  the near-tree measure — and a 4-row table (`ac`, `ac-s2m2`, `cmg-k`,
  `cmg-k-elim`) with `relres ≤ 1e-8` for all rows.
- The story to check: `cmg-k` on these graphs historically needs ~35–64 outer
  iterations; `cmg-k-elim` should cut iterations and `tot_s` sharply if the
  Spielman graphs are as tree-like as expected (small `core%`).
- Afterwards, restore any stashed local results if wanted: `git stash list`,
  `git stash pop`.

## Explicitly avoided pitfalls (why each step is shaped this way)

| Pitfall | Guard |
|---|---|
| Login OOM on precompile | No Julia on login at all; C2 runs on compute |
| Compute has no internet | `JULIA_PKG_OFFLINE=true`; demo run with `--offline`; matrices verified cached in A3 |
| Module clobbers depot (`:`) | Depot exported after `module load`, checked by echo |
| Old URL-pinned CMG silently used | C3 `pathof` + `EliminatedHierarchy` gate before any benchmarking |
| Dirty tree blocks pull | A1 stash step |
| Stale `sk<k>i<i>` names | merged FlowIPM22 loaders; A3 checks the new `Spielman_k*` names |
| k500/k600 memory blowup | medium scale for interactive; 200G batch job for paper |
