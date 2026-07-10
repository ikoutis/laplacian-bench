# laplacian-bench

A benchmark suite for **Laplacian / SDDM linear-system solvers**, comparing the
approximate-Cholesky solvers of Gao–Kyng–Spielman
([arXiv 2303.00709](https://arxiv.org/abs/2303.00709)) against
**Combinatorial Multigrid (CMG)** variants — including CMG with degree-1/2
elimination — over the paper's full problem collection. This repository is a
fork of the authors' [SDDM2023](https://rjkyng.github.io/SDDM2023) benchmark,
extended with the CMG columns, a reproducible runner, and public result tables.

## Results

Full tables: **[paper_comparison.md](performance-analyses/chol-vs-kcycle/paper_comparison.md)**
(per-instance, all solvers) · [paper_comparison.csv](performance-analyses/chol-vs-kcycle/paper_comparison.csv)
(machine-readable) · [coverage.txt](performance-analyses/chol-vs-kcycle/coverage.txt)
(completeness audit).

Five solvers over all **14 family groups (SPE included)**, every solver run on
the **same node** for each instance. Numbers are median per-family **speedups
relative to `ac`** (the paper's ApproxChol baseline, ≡ 1.00×; >1 = faster than
`ac`). **total** = build + solve:

| family | instances | `ac-s2m2` | `cmg-legacy` | `cmg-k` | `cmg-k-elim` |
|---|---:|---:|---:|---:|---:|
| sachdeva_star | 15 | 2.02× | **8.13×** | 6.80× | **6.91×** |
| checkered | 7 | 0.67× | 3.82× | **4.07×** | **3.84×** |
| aniso | 7 | 0.57× | 2.40× | 2.79× | **2.84×** |
| wgrid | 7 | 0.61× | 2.63× | **3.29×** | **2.83×** |
| uniform_grid | 3 | 0.50× | 2.21× | 2.58× | **2.65×** |
| spe | 5 | 0.69× | 1.80× | **2.51×** | **2.28×** |
| suitesparse (Laplacian) | 3 | 0.77× | 1.11× | 1.38× | **1.51×** |
| wted_chimera | 4 | 0.60× | 1.28× | 1.30× | **1.33×** |
| chimeraIPM | 128 | 0.47× | 1.23× | 1.15× | **1.20×** |
| uni_bndry_chimera | 4 | 0.66× | 1.31× | 1.32× | **1.15×** |
| uni_chimera | 4 | 0.60× | 1.16× | 1.10× | **1.10×** |
| suitesparse (SDDM) | 25 | 0.75× | 1.04× | 1.08× | **1.02×** |
| spielmanIPM | 61 | 0.65× | 0.13× | 0.20× | **0.99×** |
| wted_bndry_chimera | 4 | 0.61× | 0.94× | 1.09× | **0.90×** |

Solve-only speedup vs `ac` (excludes the one-time preconditioner build — the
number that matters when one factorization serves many right-hand sides):

| family | `ac-s2m2` | `cmg-legacy` | `cmg-k` | `cmg-k-elim` |
|---|---:|---:|---:|---:|
| sachdeva_star | 12.44× | **20.03×** | 12.84× | 12.79× |
| checkered | 0.94× | 2.72× | **3.03×** | 2.80× |
| spielmanIPM | 1.01× | 0.04× | 0.06× | **6.80×** |
| wgrid | 1.03× | 1.28× | **1.79×** | 1.59× |
| aniso | 0.96× | 1.16× | 1.39× | **1.43×** |
| uniform_grid | 1.07× | 1.00× | **1.32×** | 1.30× |
| spe | **1.19×** | 0.65× | **1.19×** | 1.04× |
| wted_chimera | **1.28×** | 0.82× | 0.92× | 1.09× |
| suitesparse (Laplacian) | **1.34×** | 0.66× | 0.97× | 0.96× |
| suitesparse (SDDM) | **1.09×** | 0.89× | 0.89× | 0.81× |
| chimeraIPM | **0.92×** | 0.49× | 0.44× | 0.58× |
| uni_chimera | **1.05×** | 0.60× | 0.64× | 0.63× |
| uni_bndry_chimera | **1.08×** | 0.59× | 0.58× | 0.58× |
| wted_bndry_chimera | **1.08×** | 0.47× | 0.62× | 0.57× |

**Summary.** On total time, `cmg-k-elim`'s per-family medians range from 6.9×
(sachdeva_star) and 2.7–3.8× (the grid families and SPE) on the structured
families down to 0.90–1.5× on the unstructured ones; its lowest family median is
`wted_bndry_chimera` at 0.90×. The columns separate the two CMG changes:

- **The K-cycle** (`cmg-k` vs `cmg-legacy`) is 10–35% faster on the grids, SPE,
  and SuiteSparse; legacy is faster on `sachdeva_star`, where the hierarchy is
  near-exact and plain PCG suffices.
- **The elimination** (`cmg-k-elim` vs `cmg-k`) is within 0.95–1.02× where there
  are no low-degree vertices to remove (the adaptive skip — grids, SPE,
  sachdeva) and larger where there are: on the near-tree `spielmanIPM` family,
  CMG without elimination runs at 0.13×/0.20× of `ac` (total), while
  `cmg-k-elim` is at parity on total time and 6.8× on solve.
- **`ac-s2m2`** is 1.4–2× slower than plain `ac` except on `sachdeva_star`
  (2.02×), where plain `ac` stagnates.
- On solve-only time `ac` and `ac-s2m2` are faster than the CMG variants on the
  unstructured families (chimeras, IPM-chimera); their totals there are
  build-dominated.

All five solvers converged on every instance. Plain `ac` reached a maximum
relres of 1.5e-6 on the largest `sachdeva_star` instances (above the 1e-8
target); the other four stayed ≤ 1e-8.

**Worst case per family** (minimum per-instance total-time speedup of
`cmg-k-elim` vs `ac` — how badly it can lose, and where):

| family | worst vs `ac` | on instance |
|---|---:|---|
| checkered | 3.26× | 2·10⁸ nnz, 32³ blocks |
| aniso | 2.28× | 2·10⁸ nnz, ξ=0.001 |
| uniform_grid | 2.27× | 2·10⁶ nnz |
| sachdeva_star | 2.05× | star_join(K₁₀₀, 50) |
| wgrid | 1.90× | 2·10⁸ nnz, w=10 |
| spe | 1.62× | spe2m |
| suitesparse (Laplacian) | 1.23× | Gaertner/nopoly |
| wted_chimera | 1.15× | n=10⁷ |
| uni_bndry_chimera | 1.07× | n≈10⁷ |
| uni_chimera | 0.98× | n=10⁶ |
| spielmanIPM | 0.78× | k600.i4 |
| chimeraIPM | 0.69× | uc.n1e5.i3.eps0.1.4 |
| wted_bndry_chimera | 0.69× | n≈10⁵ |
| suitesparse (SDDM) | 0.55× | HB/bcsstm21 |

The lowest per-instance ratios in the suite are on a small SuiteSparse mass
matrix (`bcsstm21`, millisecond scale, 0.55×) and single chimera instances
(0.69×). On the structured families and SPE the worst instance stays ≥ 1.6×.

**Solver columns.** `ac` = ApproxChol, the paper's base solver
(`ApproxCholParams(:deg,0,0)`) — the baseline (1.00×); `ac-s2m2` = its
split-2/merge-2 variant (the paper's AC-s2m2, recommended there for hard
problems); `cmg-legacy` = classic CMG (the original stationary cycle + PCG,
`cmg-v` in the raw tables); `cmg-k` = CMG K-cycle without elimination;
`cmg-k-elim` = CMG K-cycle with exact degree-1/2 elimination and the adaptive
skip — the default of
[CombinatorialMultigrid.jl](https://github.com/ikoutis/CombinatorialMultigrid.jl).

### Caveats

- **Environment.** All numbers are from compute nodes of NJIT's Wulver cluster,
  **single-threaded** (Julia 1.11.9, `JULIA_NUM_THREADS=1`, BLAS pinned to 1
  thread), tolerance 1e-8 (SPE 5e-9), maxits 1000, `b` random and projected to
  the range space — the paper's protocol. **One repetition** per instance after
  a full JIT warm-up of every solver, seed 2. Crucially, **all five solvers run
  sequentially on the same node for every instance**, so the ratios are
  node-consistent; absolute times still vary across Wulver's heterogeneous
  nodes, so don't compare seconds across families or runs. Treat single-sample
  ratios within ±15% of 1 as noise.
- **Chimera rows are one random graph per size** (the generators draw a fresh
  graph per repetition, and this run used one repetition).
- This comparison covers the two ApproxChol solvers and the CMG variants only —
  not the paper's full external-solver sweep (LAMG, HyPre, PETSc, ICC,
  MATLAB-CMG).

## The benchmark

13 problem families reproduce the paper's collection with matching generators
and sweeps: deterministic SDDM grids (`uniform_grid`, `aniso`, `wgrid`,
`checkered`), `sachdeva_star`, a curated 28-matrix SuiteSparse selection
(split Laplacian/SDDM as in the paper), the SPE reservoir-simulation matrices
(manual data), two interior-point-method families (`chimeraIPM`, `spielmanIPM`),
and four random chimera families at sizes 10⁴–10⁷. The harness records, per
instance and solver: build time, solve time, total, PCG iterations, and the
final relative residual.

## Reproducing the results (Linux)

**Requirements:** Julia ≥ 1.10, git, tens of GB of disk for the full data set.
RAM depends on scale: `smoke` and `medium` run on a laptop (≤ 32 GB);
the full `paper` scale includes instances with up to ~2.5·10⁸ nonzeros
(large grids, 10⁷-node chimeras, the biggest IPM graphs) that need
**up to ~200 GB RAM** and multi-day total compute.

```bash
# 1. Clone and set up the Julia environment (fetches CombinatorialMultigrid.jl)
git clone https://github.com/ikoutis/laplacian-bench.git
cd laplacian-bench
julia --project=. setup.jl

# 2. Download the automatic data (SuiteSparse selection + an IPM subset)
cd performance-experiments
julia --project=.. download_data.jl --scale paper     # or: medium | smoke

# 3. Manual data (required to reproduce the tables):
#    - ipmMat.zip (Tutorial.md §IPM): unzip in matrix-files/. It extracts into
#      matrix-files/ipmMat/; the loaders read it there directly.
#        cd ../matrix-files && unzip /path/to/ipmMat.zip
#    - spe.zip (Tutorial.md §SPE): spe{0.5m,2m,4m,8m,16m}.mm must sit flat in
#      matrix-files/. Without it the spe family is skipped and reported as such.
#        unzip /path/to/spe.zip && cd ../performance-experiments

# 4. Run the benchmark (sequential; the published tables used exactly this:
#    five solvers, one repetition after a full JIT warm-up, seed 2)
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./run_paper_comparison.sh run --scale paper       # or: medium | smoke

# 5. Aggregate into the public tables (same env so it picks up this run's files)
PAPER_SOLVERS="ac,ac-s2m2,cmg-v,cmg-k,cmg-k-elim" CVK_REPS=1 CVK_SEED=2 \
    ./run_paper_comparison.sh summarize
```

The outputs land in `performance-analyses/chol-vs-kcycle/`:
`paper_comparison.csv`, `paper_comparison.md`, and `coverage.txt` — the last one
audits expected-vs-produced instance counts per family so you can verify nothing
was silently skipped. A fast, data-free sanity check of the aggregation pipeline
is `./run_paper_comparison.sh selftest`.

To try a subset first, `--scale smoke` runs a small instance of every family in
minutes; individual families can be run directly, e.g.
`julia --project=.. chol_vs_kcycle.jl aniso --scale medium`.

(For cluster execution with Slurm, see
`performance-experiments/wulver/README-paper-comparison.md`.)

## Repository layout

- `julia-files/` — problem generators, loaders, and the solver harness
  (from SDDM2023, plus the CMG solver columns in `cmgSolvers.jl`).
- `performance-experiments/` — the runner (`chol_vs_kcycle.jl`,
  `run_paper_comparison.sh`), data prefetch (`download_data.jl`), and the
  table generator (`make_paper_tables.jl`).
- `performance-analyses/` — results; the public tables live in
  `chol-vs-kcycle/`.
- `matrix-files/` — downloaded/staged matrices (data files untracked; only the
  curated SuiteSparse selection list ships with the repo).

## Credits

The benchmark design, generators, and the ApproxChol solvers are due to
Yuan Gao, Rasmus Kyng, and Daniel A. Spielman —
[*Robust and Practical Solution of Laplacian Equations by Approximate
Elimination*](https://arxiv.org/abs/2303.00709); original benchmark:
[SDDM2023](https://rjkyng.github.io/SDDM2023). CMG:
[CombinatorialMultigrid.jl](https://github.com/ikoutis/CombinatorialMultigrid.jl)
(Koutis–Miller–Tolliver combinatorial multigrid, with K-cycle and degree-1/2
elimination extensions).
