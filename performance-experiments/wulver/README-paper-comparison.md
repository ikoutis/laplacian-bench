# Paper-comparison run: ApproxChol vs CMG on the full benchmark

This runbook reproduces the benchmark of **arXiv 2303.00709** (Gao–Kyng–Spielman,
*Robust and Practical Solution of Laplacian Equations by Approximate Cholesky
Factorization*) and adds the CMG variants, producing **public CSV + Markdown
tables** analogous to the paper's per-family tables.

## What runs

Four solver columns (the public table), pinned by `run_paper_comparison.sh`:

| column        | solver                                             |
|---------------|----------------------------------------------------|
| `ac`          | Laplacians.jl ApproxChol `ApproxCholParams(:deg,0,0)` — paper's default |
| `ac-s2m2`     | ApproxChol `ApproxCholParams(:deg,2,2)` — paper's split=2/merge=2 |
| `cmg-v`       | CombinatorialMultigrid legacy V-cycle / PCG        |
| `cmg-k-elim`  | CMG K-cycle **with** degree-1/2 elimination preprocessing |

`ac` and `ac-s2m2` are exactly the paper's two headline solvers. Everything else
(tol 1e-8, maxits 1000, single-threaded, `b = randn` projected to range and
mean-centered, residual `‖Ax−b‖/‖b‖`, reps = 3) matches the paper.

**Families (all 13 paper families):** `uniform_grid, aniso, wgrid, checkered,
sachdeva_star, suitesparse, chimeraIPM, spielmanIPM, uni_chimera,
uni_bndry_chimera, wted_chimera, wted_bndry_chimera` — plus **`spe` is
intentionally skipped** (see below).

## Manual data to stage first (login node)

Compute nodes are offline, so all data must be present in `matrix-files/` before
you submit. `setup_wulver.sh` auto-downloads the SuiteSparse selection and the
FlowIPM22 IPM **subset**, but two sources are manual:

1. **Full IPM sweep — REQUIRED for this comparison** (the paper's full set, not
   the FlowIPM22 subset). Download `ipmMat.zip`
   (<https://www.dropbox.com/s/qvobilehu9vzeqm/ipmMat.zip?dl=0>, Tutorial.md
   §"IPM") and unzip it. The archive contains Matrix Market `.mm` Laplacians
   named `uni_chimera.n<n>.i<i>.eps<eps>.<cnt>.mm` (chimera-IPM) and
   `spielman.k<k>.low<lo>.up<up>.i<i>.mm` (Spielman-IPM):

   ```bash
   cd /path/to/laplacian-bench/matrix-files
   unzip /path/to/ipmMat.zip                  # -> creates an ipmMat/ subdirectory
   mv ipmMat/*.mm . && rmdir ipmMat           # flatten (optional; loaders also read ipmMat/)
   ls uni_chimera.*.mm spielman.*.mm | head   # sanity check
   ```

   Note the zip extracts into an **`ipmMat/` subdirectory**. The loaders scan
   both `matrix-files/` and `matrix-files/ipmMat/`, so flattening is optional.
   (The paper's own legacy scripts expect a different short naming, `uc.i*.mm` /
   `sk*i*.mm`, from a separate SuiteSparse route — not needed here; the
   paper-comparison loaders read the Dropbox names directly.)

   When these files are present, `chimeraIPMInstances` / `spielmanIPMInstances`
   (in `julia-files/benchFamilies.jl`) automatically prefer the full sweep over
   the auto-downloaded FlowIPM22 subset. If `ipmMat.zip` is absent the loaders
   silently fall back to the subset, so **stage it or the IPM tables will be the
   thin FlowIPM22 subset instead of the paper's full sweep.**

2. **SPE — SKIPPED** (user decision). The SPE matrices (`spe.zip`,
   ~0.5m–16m variables) are not staged, so the `spe` family produces no rows.
   `coverage.txt` lists it as `skipped` rather than missing. To include it later,
   download `spe.zip` (Tutorial.md §"SPE benchmark") into `matrix-files/` and add
   `spe` back — nothing else changes.

## Procedure

All commands are run from `performance-experiments/` unless noted.

### 1. Setup (login node — has internet)

```bash
cd performance-experiments/wulver
# stage ipmMat.zip into ../../matrix-files first (see above)
./setup_wulver.sh --scale paper        # Julia env + SuiteSparse + IPM subset prefetch
```

If your `$HOME` quota is tight, `export CVK_DEPOT=/project/ikoutis/$USER/.julia`
(or another roomy path) before running — the sbatch scripts honor the same var.

**Interactive Julia (login or compute node).** To run Julia by hand — e.g. to
confirm the IPM sweep is seen — source the env helper instead of retyping the
module/depot loads (compute nodes are offline, so `srun` onto one first):

```bash
srun --account=ikoutis --partition=general --qos=standard \
     --cpus-per-task=4 --mem=16G --time=00:30:00 --pty bash
source performance-experiments/wulver/env_wulver.sh     # module + depot + thread pins
julia --project=. -e 'include("julia-files/benchFamilies.jl"); @show length(chimeraIPMInstances(:paper)), length(spielmanIPMInstances(:paper))'
```

(Run `setup_wulver.sh` once on the login node first so the depot is instantiated;
the helper only sets the environment, it does not install packages.)

### 2. Submit the full run (login node → Slurm)

```bash
cd performance-experiments
CVK_REPS=3 ./run_paper_comparison.sh submit --scale paper --account ikoutis
```

This pins `--solvers ac,ac-s2m2,cmg-v,cmg-k-elim` and reps=3, emits the
small/large task manifests, and submits two array jobs (the large tier gets
`--mem=200G --time=72:00:00` for the big grids, 1e7 chimeras, and k500/k600
Spielman-IPM graphs). Monitor with `squeue -u $USER`; per-family `.jld2` results
land in `../performance-analyses/chol-vs-kcycle/`.

**Single-allocation alternative** (smoke/medium, or a modest paper run on one
big-memory node): edit `--account` in `paper_comparison.sbatch`, then

```bash
sbatch --export=ALL,CVK_SCALE=medium paper_comparison.sbatch
# or, paper scale on one big node:
sbatch --mem=200G --time=72:00:00 --export=ALL,CVK_SCALE=paper paper_comparison.sbatch
```

It runs the whole suite sequentially and summarizes at the end.

### 3. Summarize into public tables (after all jobs finish)

```bash
cd performance-experiments
./run_paper_comparison.sh summarize
```

Writes into `../performance-analyses/chol-vs-kcycle/`:

- **`paper_comparison.csv`** — machine-readable, RFC-4180 quoted (some testNames
  contain commas): one row per instance (all instances, unfiltered), the system
  `kind` (lap/sddm), and per solver `build/solve/tot/its/err` medians + failure
  count.
- **`paper_comparison.md`** — one table per family, median total seconds with
  median iterations in parentheses, a per-family max-relres line, and a median
  `cmg-k-elim` vs `ac` total-time speedup footer. SuiteSparse is split into a
  Laplacian and an SDDM table with `ne > 1000`, mirroring the paper. Analogous
  to the paper's tables.
- **`coverage.txt`** — expected vs produced instance counts per family and
  per-solver failures, so **"nothing was missed" is auditable** (IPM full sweep
  present; SPE listed as intentionally skipped).

The `.csv` and `.md` are the artifacts to share publicly.

## Verifying nothing was missed

Open `coverage.txt`. Each fixed-sweep family should read `OK` at its paper count
(`uniform_grid 3, aniso 7, wgrid 7, checkered 7, sachdeva_star 15,
suitesparse 28`, chimera families 4 each). `chimeraIPM`/`spielmanIPM` report the
produced sweep count (data-dependent on `ipmMat.zip`); `spe` reads `skipped`.
Any `PARTIAL`/`MISSING` flags a family to re-run or stage data for. The
per-solver failure lines flag instances where a solver hit maxits or errored.

## Local gates (no cluster, no data)

```bash
# summarizer correctness (CSV quoting + median aggregation), no data needed:
./run_paper_comparison.sh selftest

# tiny end-to-end shape check on real (smoke) output:
./run_paper_comparison.sh run --scale smoke
./run_paper_comparison.sh summarize
```
