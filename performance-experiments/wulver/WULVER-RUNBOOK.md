# Wulver runbook — the workflow we follow

Consolidated, session-to-session steps for running this benchmark (and a modified
`CombinatorialMultigrid.jl`) on NJIT's **Wulver** cluster. Two independent
concerns are separated below: **Part A — System** (Slurm / login-vs-compute /
how to get work onto a node) and **Part B — Julia** (env, depot, the dev-linked
CMG, running things). Part C is the concrete recurring loop; Part D is the
templates for sending jobs to a compute node **without** an interactive "grab".

Task-specific runbooks live alongside this file and are referenced where needed:
`setup_wulver.sh` (one-time setup), `README-paper-comparison.md` (the full paper
run), `UPDATING-CMG-ON-WULVER.md` (running a modified CMG), `env_wulver.sh` (the
source-able environment).

Account/paths used throughout (edit for another user):
`--account=ikoutis`, repo at `/project/ikoutis/github/laplacian-bench`, Julia
depot `/project/ikoutis/$USER/.julia`.

---

## Part A — System (Slurm / Wulver)

### A.1 The one rule: login vs compute
- **Login node** — has **internet** and is where you do `git`, edit files, and
  submit jobs. Its memory cgroup is **too tight to precompile Julia** (it
  OOM-kills — even a first-time `using JLD2`). Never run Julia work of substance
  here.
- **Compute node** — has real RAM **and** outbound internet. Do **all** Julia
  work here: precompile, `Pkg` operations, running scripts, `Pkg.test`,
  summarizing.

So every session is: **git/setup on login → run on compute → commit on login.**

### A.2 Getting work onto a compute node — four ways
| Way | Command shape | Use when |
|---|---|---|
| Interactive | `srun … --pty bash` | exploring / debugging, iterating by hand |
| One-shot foreground | `srun … bash -lc '…'` | a known command, want live output, will wait |
| Detached batch | `sbatch … --wrap='…'` | **preferred for real runs** — logged, survives logout |
| Array job | `sbatch --array=… chol_vs_kcycle_array.sbatch` | the benchmark sweep (many tasks) |

Templates for the non-interactive ways are in **Part D** — prefer them over
"grab a node then type"; the interactive `--pty` is only for exploration.

### A.3 Resource tiers (what to pass)
- **Small** (most jobs): `--mem=16G` (up to `48G`) `--cpus-per-task=4`
  `--time=01:00:00`…`24:00:00`.
- **Large** (2e8-nnz grids, 1e6/1e7 chimeras, spielmanIPM k500/k600):
  `--mem=200G --time=48:00:00` (bump to `72:00:00` for the biggest un-chunked
  graphs).
- Solvers are single-threaded (see B.1), so `--cpus-per-task=4` is plenty; the
  extra cores just cap parallel precompile.

### A.4 Common Slurm flags (fixed for this cluster)
```
--account=ikoutis --partition=general --qos=standard --nodes=1 --ntasks=1
```

### A.5 Monitor / manage
```bash
squeue -u $USER                       # queued + running jobs
scancel <jobid>                       # cancel
ls -lt performance-experiments/wulver/logs/   # sbatch output (%A_%a.out)
```

### A.6 Git (login node only — needs internet; SSH key has a passphrase)
```bash
cd /project/ikoutis/github/laplacian-bench
git fetch origin <branch> && git checkout <branch> && git pull origin <branch>
```
Results are committed/pushed from the login node after a run (Part C step 4).

---

## Part B — Julia

### B.1 The environment: source `env_wulver.sh`
```bash
source performance-experiments/wulver/env_wulver.sh
```
It sets, matching the sbatch scripts:
- **module + Julia**: `module load Julia/1.11.9` (NOT the bare default 1.12.6;
  there is no 1.10 module).
- **depot**: `JULIA_DEPOT_PATH=/project/ikoutis/$USER/.julia` (keeps packages off
  the small `$HOME` quota; must match what `setup_wulver.sh` used).
- **offline**: `JULIA_PKG_OFFLINE=true` (compute nodes shouldn't reach the
  registry mid-run; override `JULIA_PKG_OFFLINE=false` on the **login** node for
  installs).
- **threads**: `JULIA_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`,
  `OMP_NUM_THREADS=1` — the paper protocol; keeps timings fair.

Source it in every compute-node shell (and on the login node before a `Pkg`
install). It only sets env — it does not instantiate/precompile.

### B.2 Which project
- From the repo root: `julia --project=.`
- From `performance-experiments/`: `julia --project=..`

### B.3 The dev-linked CombinatorialMultigrid.jl (important)
The benchmark does **not** use a registered CMG — `setup.jl` `Pkg.develop`s a
**local checkout** (via `CMG_DEV`), so the environment loads CMG from a path on
disk. In our sessions that checkout is **nested inside the repo** at
`/project/ikoutis/github/laplacian-bench/CombinatorialMultigrid.jl`.

Consequences you must remember:
- **To run updated CMG, pull *that* checkout** — not some other clone:
  ```bash
  cd /project/ikoutis/github/laplacian-bench/CombinatorialMultigrid.jl
  git fetch origin <branch> && git checkout <branch> && git pull origin <branch>
  ```
  Because it's a `Pkg.develop` path dependency, Julia **recompiles it from the
  updated source on the next `using`** — no `setup.jl` re-run needed for source
  edits.
- **Verify the right code is loaded** (the fastest sanity check):
  ```bash
  julia --project=. -e 'using CombinatorialMultigrid; println(pathof(CombinatorialMultigrid))'
  ```
  It should print the `…/laplacian-bench/CombinatorialMultigrid.jl/src/…` path,
  not a `$JULIA_DEPOT_PATH/packages/CombinatorialMultigrid/<hash>/…` one.
- First-time / from-scratch selection of a CMG version is done via `setup.jl`
  with `CMG_DEV=<path>` (local checkout) or `CMG_REV=<commit-or-branch>` (pushed
  git rev). See `UPDATING-CMG-ON-WULVER.md`.

### B.4 Resolve / instantiate / precompile
- **`Pkg.resolve()` / `Pkg.instantiate()`** touch the registry, so run them on the
  **login** node when a `Project.toml` changed (e.g. a dep was added/removed):
  ```bash
  # login node, once:
  julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
  ```
- **Precompilation is heavy and memory-hungry** → it happens on first `using` on
  the **compute** node (JLD2 especially). Let it run; it caches per-CPU.

### B.5 Running things (compute node)
```bash
source performance-experiments/wulver/env_wulver.sh
cd performance-experiments
julia --project=.. <script>.jl <args>          # e.g. compare_split.jl, diagnose_*.jl
# CMG unit tests:
julia --project=/project/ikoutis/github/laplacian-bench/CombinatorialMultigrid.jl -e 'using Pkg; Pkg.test()'
```

---

## Part C — The recurring loop (what we actually do each session)

1. **Login — update code.**
   ```bash
   cd /project/ikoutis/github/laplacian-bench
   git fetch origin <branch> && git checkout <branch> && git pull origin <branch>
   # if you changed which CMG to run, also pull the dev checkout (B.3):
   cd CombinatorialMultigrid.jl && git checkout <branch> && git pull && cd ..
   # only if a Project.toml changed:
   julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
   ```
2. **Compute — run** (interactively per A.2, or submit per Part D). Always
   `source env_wulver.sh` first.
3. **Compute — summarize** (if applicable):
   ```bash
   PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
     ./run_paper_comparison.sh summarize
   ```
4. **Login — commit results.**
   ```bash
   git add performance-analyses/chol-vs-kcycle/paper_comparison.{csv,md} \
           performance-analyses/chol-vs-kcycle/coverage.txt
   git commit -m "…" && git push origin <branch>
   ```

---

## Part D — Sending a job to a compute node without "grabbing" one

Prereq (once, when deps changed): `Pkg.resolve(); Pkg.instantiate()` on the
**login** node (D-jobs run offline). The env is sourced **inside** the job so it
gets the module/depot on the compute node.

### D.1 Detached batch (recommended) — `sbatch --wrap`
Runs on a compute node, logs to a file, survives logout. No separate `.sbatch`.
```bash
cd /project/ikoutis/github/laplacian-bench/performance-experiments
mkdir -p wulver/logs
sbatch --account=ikoutis --partition=general --qos=standard \
       --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=01:00:00 \
       --job-name=cvk_oneoff --output=wulver/logs/cvk_oneoff_%j.out \
       --wrap='source wulver/env_wulver.sh && julia --project=.. compare_split.jl'
squeue -u $USER
tail -f wulver/logs/cvk_oneoff_*.out
```
Bump `--mem`/`--time` (and add nothing else) for a bigger job.

### D.2 One-shot foreground — `srun` with a command (no `--pty`)
Allocates a node, runs the command there with **live output**, releases it. Blocks
your shell until done.
```bash
cd /project/ikoutis/github/laplacian-bench/performance-experiments
srun --account=ikoutis --partition=general --qos=standard \
     --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=16G --time=01:00:00 \
     bash -lc 'source wulver/env_wulver.sh && julia --project=.. compare_split.jl'
```

### D.3 The benchmark sweep — array job
For the paper run / large sweeps, don't hand-roll; use the existing flow, which
submits `chol_vs_kcycle_array.sbatch` for you:
```bash
CVK_REPS=3 ./run_paper_comparison.sh submit --scale paper --account ikoutis
```
(small tier + a `--mem=200G --time=48:00:00` large tier; see
`README-paper-comparison.md`).

### D.4 When to still use interactive `srun … --pty bash`
Only for exploration/debugging where you don't yet know the command — poking at a
graph, iterating on a script, reading errors. For anything you can name up front,
prefer D.1 (detached, logged) or D.2 (foreground).

---

## Part E — Sparsify-on-stall benchmark comparison

Compares `ac`, `cmg-k-elim` (reference), and the two sparsify-on-stall columns
`cmg-sparsify-l` (:legacy) and `cmg-sparsify-ks` (:kscycle), which inject a
spanner + uniform spectral-sparsifier level when CMG aggregation stalls (both
build with `eliminate=true`). New `_inj` column records injected levels per
solver per instance; `make_paper_tables.jl` emits a **Stall-triggered
comparison** section listing every instance where sparsify fired (`inj > 0`).

### E.1 Pin the CMG branch (login node)
The benchmark's dev-linked CMG (B.3) must be on the feature branch, which adds a
`DataStructures` dependency — so re-resolve:
```bash
cd /project/ikoutis/github/laplacian-bench/CombinatorialMultigrid.jl
git fetch origin claude/julia-sparsify-on-stall
git checkout claude/julia-sparsify-on-stall && git pull origin claude/julia-sparsify-on-stall
cd /project/ikoutis/github/laplacian-bench
git fetch origin claude/cmg-degree-elimination-qwi3iq
git checkout claude/cmg-degree-elimination-qwi3iq && git pull origin claude/cmg-degree-elimination-qwi3iq
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'   # picks up DataStructures
```
Then precompile once on a **compute** node (B.4) and confirm `pathof` points at
the nested checkout.

The four columns for every run:
```bash
export PAPER_SOLVERS="ac,cmg-k-elim,cmg-sparsify-l,cmg-sparsify-ks"
```

### E.2 Run 1 — chimeras at 1e6 (the known-stall set)
```bash
cd /project/ikoutis/github/laplacian-bench/performance-experiments
CVK_ONLY="uni_chimera uni_bndry_chimera wted_chimera wted_bndry_chimera" \
CVK_CHIMERA_SIZES="1e6" \
PAPER_SOLVERS="$PAPER_SOLVERS" CVK_REPS=3 \
  ./run_paper_comparison.sh submit --scale paper --account ikoutis
```
1e6 chimera chunks route to the large tier (`--mem=200G --time=48:00:00`,
`--chunk` core-parallel, whole instances per node → same-node-per-instance).

### E.3 Run 2 — other families, stall-triggered only
Run the remaining families; the comparison is only meaningful where sparsify
fired, which the Stall-triggered table surfaces (rows with `inj = 0` stay in the
full per-family tables but are not in that section):
```bash
CVK_ONLY="uniform_grid aniso wgrid checkered sachdeva_star suitesparse spe chimeraIPM spielmanIPM" \
PAPER_SOLVERS="$PAPER_SOLVERS" CVK_REPS=3 \
  ./run_paper_comparison.sh submit --scale paper --account ikoutis
```

### E.4 Run 3 — artificial dense instances (deterministic stalls)
```bash
CVK_ONLY="dense_blob" PAPER_SOLVERS="$PAPER_SOLVERS" CVK_REPS=3 \
  ./run_paper_comparison.sh submit --scale paper --account ikoutis
```
`dense_blob` is small (≤1500-node blobs / blob chains from the CMG-python port);
it stays on the small tier and every instance stalls, so all four columns run.

### E.5 Summarize + commit
```bash
# compute node, after the jobs finish:
julia --project=.. make_paper_tables.jl \
  --solvers ac,cmg-k-elim,cmg-sparsify-l,cmg-sparsify-ks
# sanity self-test of the table/inj/stall logic (no data needed):
julia --project=.. make_paper_tables.jl --selftest
```
Then commit `performance-analyses/chol-vs-kcycle/{paper_comparison.csv,paper_comparison.md,coverage.txt}`
per Part C step 4. Read the **Stall-triggered comparison** section and the mean
injected-levels line to see where and how hard sparsify fired.

---

## Gotchas (the ones that bit us)
- **Precompile only on compute** — login OOM-kills it.
- **`module load Julia/1.11.9` explicitly** — bare `Julia` is 1.12.6.
- **`JULIA_DEPOT_PATH=/project/ikoutis/$USER/.julia`** must be set on login and in
  jobs (env_wulver.sh / the sbatch scripts do it) or Julia won't find packages.
- **Pull the *nested* CMG checkout** (`laplacian-bench/CombinatorialMultigrid.jl`)
  to change what CMG runs; confirm with `pathof(CombinatorialMultigrid)`.
- **Chunked result filenames** carry `.chunkKofC` before `.seedS.repsR.jld2`; a
  targeted re-run must remove the stale ones for that size first (see the top-up
  section of `README-paper-comparison.md`).
