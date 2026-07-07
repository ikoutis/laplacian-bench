# Running an updated CombinatorialMultigrid.jl on Wulver

Handoff notes for building and running a **modified** `CombinatorialMultigrid.jl`
(e.g. one that adds degree-1/degree-2 node elimination) against this benchmark
on NJIT's **Wulver** cluster. Written from a session that got the stock CMG
running end-to-end; the environment gotchas below are the ones that actually
bit us.

## 0. Why this matters (the motivation for the CMG change)

On the FlowIPM22 **Spielman** graphs, stock CMG converges but is slow (35–64
iterations) while approximate-Cholesky (AC) needs only ~7. Reason: those graphs
are **near-trees** — a spanning tree plus only `O(k²)` off-tree edges (off-tree
fraction ≈ `1.5/k` of the vertices, i.e. 1.5% at k100 → 0.25% at k600). A
tree-preserving elimination (AC) is near-optimal there; aggregation multigrid is
not. **Eliminating degree-1 (leaves) and degree-2 (chain) nodes collapses
exactly the tree/path parts**, so the updated CMG should slash iterations/build
on Spielman/near-tree instances. That's the change to validate here.

## 1. Where CMG lives and how the benchmark selects it

- The benchmark calls CMG through `julia-files/cmgSolvers.jl` →
  `cmg_preconditioner_lap(sys; …)` (build) and `cmg_solve(H, b; …)` (solve).
  `make_cmg_lap` builds on `lap(a)`; `make_cmg_sddm` passes the SDDM matrix and
  CMG grounds strictly-dominant rows with one extra coordinate.
- The **environment** pins CMG in `setup.jl` (repo root):
  ```julia
  const CMG_DEFAULT_REV = "08c515ed76aabc4caef08abc8bc60b98a549ea1e"
  if haskey(ENV, "CMG_DEV")
      Pkg.develop(path = ENV["CMG_DEV"])          # local checkout — picks up edits
  else
      rev = get(ENV, "CMG_REV", CMG_DEFAULT_REV)
      Pkg.add(url = "https://github.com/ikoutis/CombinatorialMultigrid.jl", rev = rev)
  end
  ```
  So the two knobs to run a modified CMG are **`CMG_DEV`** (a local source
  checkout via `Pkg.develop` — best for active development, edits are picked up
  on the next precompile) and **`CMG_REV`** (a pushed git branch/commit).
- Stock installed source on the Wulver depot lives at
  `$JULIA_DEPOT_PATH/packages/CombinatorialMultigrid/<hash>/src/` (main algorithm
  `cmgAlg.jl`; `cmg_preconditioner_lap` / `cmg_solve` are the public entry
  points). A `Pkg.develop`'d version instead loads from your checkout path.

## 2. Wulver environment (always set these)

```bash
module load Julia/1.11.9          # NOT the default 1.12.6, and there is no 1.10 module
export JULIA_DEPOT_PATH=/project/ikoutis/$USER/.julia   # keep the depot off the ~50GB $HOME quota
```
Repo: `/project/ikoutis/github/laplacian-bench`
(branch `claude/julia-wulver-hpc-benchmark-m2qh8r`).

## 3. LOGIN vs COMPUTE — the #1 gotcha

- **Login node**: has internet, **but its memory cgroup is too tight to
  precompile Julia** — precompilation gets OOM-killed (`Killed`). Do **not** run
  `setup.jl`, `Pkg.precompile`, or any heavy `using` on the login node. Even
  reading a `.jld2` via Julia can trip it (first-time `JLD2` precompile).
- **Compute nodes**: have real RAM (precompile succeeds) **and** outbound
  internet (confirmed — they can `Pkg.add` from GitHub and download matrices).
  Precompile caches are **CPU-specific**, so a fresh compute node may rebuild
  ~35 packages the first time — that's fine, it has the memory.
- **Rule: do all serious work — building the modified CMG, precompiling, and
  running it — on a COMPUTE node.** Use the login node only for `git`, editing
  files, and `sbatch`/`squeue`.

## 4. Build + run the updated CMG (interactive, on a compute node)

```bash
# 1. get a compute node (RAM + internet):
srun --account=ikoutis --partition=general --qos=standard \
     --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=48G --time=02:00:00 --pty bash

# 2. environment:
module load Julia/1.11.9
export JULIA_DEPOT_PATH=/project/ikoutis/$USER/.julia
export JULIA_NUM_PRECOMPILE_TASKS=8     # cap parallel precompile workers if memory is tight
cd /project/ikoutis/github/laplacian-bench

# 3a. OPTION A — local dev checkout (recommended while iterating on the source):
#     put your modified CombinatorialMultigrid.jl somewhere on /project, e.g.:
#       git clone https://github.com/ikoutis/CombinatorialMultigrid.jl \
#                 /project/ikoutis/github/CombinatorialMultigrid.jl
#       cd /project/ikoutis/github/CombinatorialMultigrid.jl && git checkout <your-branch>
#     then point the benchmark env at it and rebuild:
CMG_DEV=/project/ikoutis/github/CombinatorialMultigrid.jl julia --project=. setup.jl

# 3b. OPTION B — a pushed branch/commit instead of a local checkout:
CMG_REV=<commit-or-branch> julia --project=. setup.jl

# 4. verify you're running the UPDATED code, not the pinned package:
julia --project=. -e 'using CombinatorialMultigrid; println(pathof(CombinatorialMultigrid))'
#     Option A should print your /project/.../CombinatorialMultigrid.jl/src/... path,
#     NOT a $JULIA_DEPOT_PATH/packages/CombinatorialMultigrid/<hash>/ path.

# 5. quick functional run — the Spielman case the change targets (k100 only, fast):
cd performance-experiments
julia --project=.. chol_vs_kcycle.jl spielmanIPM --scale smoke --reps 1 --seed 1
#     watch the cmg-k / cmg-v `iter:` and `err:` — the change should REDUCE the
#     iteration counts on Spielman vs stock CMG (stock was 35 at k100), err ≤ 1e-8.
```

### Iterating on edits (Option A)
`Pkg.develop` means Julia recompiles the package from your checkout whenever its
source changes. After editing the CMG source you do **not** re-run `setup.jl`;
just re-precompile on a compute node and run:
```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

### Reverting to stock CMG
```bash
julia --project=. -e 'using Pkg; Pkg.free("CombinatorialMultigrid")'   # undo develop
julia --project=. setup.jl                                             # back to CMG_DEFAULT_REV
```

## 5. Running at scale — use sbatch, not srun

Interactive `srun` is for quick checks. For the full Spielman sweep (k500/k600
are 126M/217M nnz) you need the **big-memory tier**; submit an array job:
```bash
cd performance-experiments/wulver
echo "spielmanIPM --scale paper --reps 1 --seed 1" > manifest.ipm.txt
sbatch --account=ikoutis --array=1-1 --mem=200G --time=72:00:00 \
       --export=ALL,CVK_MANIFEST=$PWD/manifest.ipm.txt chol_vs_kcycle_array.sbatch
```
The sbatch already `module load`s `Julia/1.11.9` and sets `JULIA_DEPOT_PATH`
(via `CVK_DEPOT`, default `/project/ikoutis/$USER/.julia`). **Batch jobs read the
project `Manifest.toml`**, so as long as you ran the `CMG_DEV=… setup.jl` step
(which rewrites the Manifest to point at your checkout), the compute-node jobs
pick up the modified CMG automatically — `/project` is shared across nodes, so
keep the checkout in place. Results stream to
`performance-analyses/chol-vs-kcycle/spielmanIPM.paper.seed1.reps1.jld2`, saved
after every rep.

## 6. Reading results without Julia-on-login

Login can't run the summarizer (precompile OOM). Either read results on a compute
node, or note the `.jld2` files are HDF5 under the hood — but they store the whole
result Dict as one serialized `dic` object (HDF5 references), so `h5py` can't
easily flatten it. Simplest: on a **compute** node,
```bash
cd /project/ikoutis/github/laplacian-bench
julia --project=. performance-experiments/summarize_chol_vs_kcycle.jl \
      performance-analyses/chol-vs-kcycle | grep -A20 spielmanIPM
```
For exact residuals, load the file and print `d["cmg-k_err"]` / `d["cmg-v_err"]`
(the dict is under key `"dic"` → `get(load(f), "dic", load(f))`).

## 7. Cheat-sheet of gotchas

- Never precompile on the login node — it OOM-kills. Compute only.
- `module load Julia/1.11.9` explicitly (bare `Julia` gives 1.12.6).
- Always `export JULIA_DEPOT_PATH=/project/ikoutis/$USER/.julia` (login + jobs).
- `pathof(CombinatorialMultigrid)` is the fastest check that the *modified* code
  is loaded.
- `Manifest.toml` after `CMG_DEV` points at a local path — don't `Pkg.instantiate`
  from a clean Manifest expecting upstream CMG; and keep the checkout on disk.
- The warning `CMG convergence may be slow due to matrix density`
  (`cmgAlg.jl:206`) is benign — it printed on some instances that still converged.
- Compute nodes have internet here, so downloads/`Pkg.add` work on them too;
  you don't have to prefetch on login if you build on compute.
