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
`--mem=200G --time=48:00:00`, which fits the chunked 10⁶/10⁷ chimeras). At paper
scale the chimera families are split into parallel `--chunk` array elements —
see *Chimera top-up* below for the chunk defaults and how whole instances stay
on one node. If a from-scratch run's large tier also carries the biggest
un-chunked fixed graphs (nnz~2e8 grids, spielmanIPM k500/k600), give those more
head-room with `--sbatch-extra "--time=72:00:00"`. Monitor with
`squeue -u $USER`; per-family `.jld2` results land in
`../performance-analyses/chol-vs-kcycle/`.

**Single-allocation alternative** (smoke/medium, or a modest paper run on one
big-memory node): edit `--account` in `paper_comparison.sbatch`, then

```bash
sbatch --export=ALL,CVK_SCALE=medium paper_comparison.sbatch
# or, paper scale on one big node:
sbatch --mem=200G --time=72:00:00 --export=ALL,CVK_SCALE=paper paper_comparison.sbatch
```

It runs the whole suite sequentially and summarizes at the end.

### 3. Summarize into public tables (after all jobs finish)

`make_paper_tables.jl` does `using JLD2, …`, and JLD2's first-use precompile is
heavy enough that the **login node may kill it** — run `summarize` on a **compute
node** (small allocation is plenty; it just reads the `.jld2` results):

```bash
srun --account=ikoutis --partition=general --qos=standard \
     --cpus-per-task=4 --mem=16G --time=00:30:00 --pty bash
source performance-experiments/wulver/env_wulver.sh
./performance-experiments/run_paper_comparison.sh summarize
```

`summarize` reads **only this run's files** — those matching
`*.seed${SEED}.reps${REPS}.jld2` in the results directory — so stale results from
earlier runs (older smoke/`reps1` dev runs, or sample `.jld2` that ship in the
repo) can't inflate the coverage counts. To summarize a different set, pass files
or a directory explicitly: `./run_paper_comparison.sh summarize <files-or-dir>`.

Writes into `performance-analyses/chol-vs-kcycle/`:

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

## Round 2: five solvers, SPE included, adaptive elimination

A second full run designed to decompose the round-1 findings:

- **Five solver columns** — `ac, ac-s2m2, cmg-v, cmg-k, cmg-k-elim`. Adding
  plain `cmg-k` separates the two effects that round 1 conflated: `cmg-k` vs
  `cmg-v` isolates the K-cycle, `cmg-k-elim` vs `cmg-k` isolates the
  elimination. The speedup summary reports both automatically.
- **Adaptive elimination** — requires the updated CombinatorialMultigrid.jl
  (`build_eliminated_hierarchy` pre-scan skip), so `cmg-k-elim`'s build no
  longer pays the rebuild on structure-less graphs. Update the in-repo CMG
  checkout before submitting.
- **reps = 1** — the harness warms up every solver on a tiny Laplacian + SDDM
  before any timed region, so JIT never lands in a timing; a single repetition
  suffices for this comparison. (Chimera families then have one random sample
  per size.)
- **seed = 2** — gives this run its own `*.seed2.reps1.jld2` filenames, so it
  cannot collide with round 1 (`seed1.reps3`) or older dev files, and
  `summarize` picks up exactly this run.
- **Same node per instance** — guaranteed by construction: each family (or
  chimera family × size) is one Slurm array task on one node, and the harness
  runs every solver sequentially on each instance within that task. No
  instance is ever timed across different nodes.
- **SPE included** — stage it first (manual): download `spe.zip`
  (<https://www.dropbox.com/s/7fp4yq69brcew8g/spe.zip?dl=0>, Tutorial.md
  §"SPE benchmark") and put `spe{0.5m,2m,4m,8m,16m}.mm` **flat** in
  `matrix-files/` (if the zip extracts into a subdirectory, `mv` them up —
  the SPE loader reads the top level only). The family's 5e-9 tolerance is
  applied automatically.

```bash
# --- login node ---
cd /project/ikoutis/github/laplacian-bench
git pull
performance-experiments/wulver/fetch_cmg.sh          # update in-repo CMG (adaptive skip)
ls matrix-files/spe*.mm                              # SPE staged? expect 5 files
cd performance-experiments
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./run_paper_comparison.sh submit --scale paper --account ikoutis

# --- after the jobs finish: compute node ---
source performance-experiments/wulver/env_wulver.sh
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./performance-experiments/run_paper_comparison.sh summarize
```

The first compute-node Julia use after `fetch_cmg.sh` precompiles the updated
CMG (offline is fine — the checkout is in-repo). Expected coverage: as round 1
plus `spe 5 present`.

## Chimera top-up: the paper's per-size instance counts

By default the chimera families now run the paper's per-size counts (seed
indices `i = 1..C`; `chimeraAndIPM.tex`): **103 / 105 / 23 / 8** instances at
sizes 10⁴ / 10⁵ / 10⁶ / 10⁷ — versus the single instance per size used in the
first two runs. The `make_paper_tables` collapse already medians per size, so
the chimera rows become per-size medians over those C samples automatically.

**Parallelism.** At `--scale paper` each chimera family × size is now split into
several Slurm array elements via `--chunk K/C` (defaults: 10⁴ → 1, 10⁵ → 2,
10⁶ → 8, 10⁷ → 8; override a size with `CVK_CHUNKS_1e6=12` etc., or all sizes
with `CVK_CHUNKS`). Each chunk runs **whole instances** — all solvers for a
given graph stay together on one node — so same-node-per-instance timing is
preserved while Slurm spreads the chunks across cores/nodes. The 10⁶/10⁷ chunks
go to the 200 G large tier at `--time=48:00:00`; the 10⁴/10⁵ chunks stay on the
small tier.

Because a chunk's filename carries a `.chunkKofC` tag
(`uni_chimera.n1e6.chunk3of8.seed2.reps1.jld2`), the new chunked files do **not**
overwrite the prior run's single-sample chimera files
(`uni_chimera.n1e6.seed2.reps1.jld2`). **Delete the stale chimera outputs first**
so the summarize glob doesn't fold an old single sample into the new per-size
median:

```bash
# --- login node ---
cd /project/ikoutis/github/laplacian-bench/performance-experiments

# 1. Clear the previous chimera results for this seed/reps (chunked or not).
#    Grids / IPM / SPE / SuiteSparse files are left untouched.
rm -f ../performance-analyses/chol-vs-kcycle/{uni_chimera,uni_bndry_chimera,wted_chimera,wted_bndry_chimera}.*.seed2.reps1.jld2

# 2. Submit the chimera top-up (chunked + 48h flow through automatically).
CVK_ONLY="uni_chimera uni_bndry_chimera wted_chimera wted_bndry_chimera" \
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./run_paper_comparison.sh submit --scale paper --account ikoutis

# --- after the chimera jobs finish: compute node, re-summarize everything ---
source performance-experiments/wulver/env_wulver.sh
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./performance-experiments/run_paper_comparison.sh summarize
```

Re-summarizing picks up the refreshed chimera chunk files together with the
untouched grid/IPM/SPE/SuiteSparse results, so only the chimera rows change.
`make_paper_tables` collapses by `(family, nv)`, so the several chunk files for
one size are medianed together into a single per-size row automatically.

**Compute note.** This is a large jump: 103 / 105 / 23 / 8 instances per family
(× 4 families × 5 solvers) instead of 4. Chunking is what keeps it tractable —
the 8× 10⁷ per family run one instance per array element, and the 23× 10⁶ run
~3 per element, all in parallel rather than as one long sequential task. Cap
counts with `CVK_EXTRA="--limit N"` per task for a faster partial top-up first,
or raise `CVK_CHUNKS_1e6`/`CVK_CHUNKS_1e7` to fan out further.
