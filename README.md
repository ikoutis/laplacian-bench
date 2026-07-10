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

Median **speedup of `cmg-k-elim`** over `ac` and over classic CMG (`cmg-v`),
per family (>1 = `cmg-k-elim` faster). **total** = build + solve; **solve**
excludes the one-time preconditioner build — the number that matters when one
factorization is reused across many right-hand sides:

| family | instances | vs `ac` total | vs `ac` solve | vs `cmg-v` total | vs `cmg-v` solve |
|---|---:|---:|---:|---:|---:|
| sachdeva_star | 15 | **4.35×** | **8.73×** | 0.72× | 0.60× |
| checkered | 7 | **2.78×** | **2.67×** | 0.86× | 1.15× |
| aniso | 7 | **2.24×** | **1.46×** | 0.94× | 1.27× |
| uniform_grid | 3 | **1.41×** | 0.76× | 0.95× | 1.27× |
| wgrid | 7 | **1.39×** | 0.83× | 0.93× | 1.21× |
| suitesparse (Laplacian) | 3 | 0.98× | 0.80× | 1.01× | 1.19× |
| spielmanIPM | 61 | 0.97× | **6.52×** | **7.98×** | **202×** |
| wted_chimera | 4 | 0.92× | 1.03× | 0.99× | 1.39× |
| suitesparse (SDDM) | 25 | 0.90× | 0.79× | 0.78× | 0.90× |
| chimeraIPM | 128 | 0.90× | 0.56× | 0.96× | 1.05× |
| uni_chimera | 4 | 0.84× | 0.80× | 1.08× | 1.43× |
| uni_bndry_chimera | 4 | 0.78× | 0.71× | 1.09× | 1.42× |
| wted_bndry_chimera | 4 | 0.76× | 0.68× | 0.93× | 1.33× |

**Summary.** On total time, CMG (K-cycle with degree-1/2 elimination) is
**1.4–4.4× faster than ApproxChol on the structured/geometric families** and
within **2–24% on the unstructured ones**. Against classic CMG it is a
**~8× total-time win on the near-tree `spielmanIPM` family** — exactly the
structure degree-1/2 elimination targets — and roughly at parity elsewhere; on
families with *no* low-degree structure classic CMG is somewhat faster (e.g.
`sachdeva_star` 0.72×), which is why elimination is a switchable option. The
solve-only columns show where the time goes: `cmg-k-elim`'s totals are
build-dominated, so in the many-right-hand-sides regime its lead grows —
**8.7× over `ac` on `sachdeva_star`, 6.5× on `spielmanIPM`** (where the
eliminated solve is ~200× faster than classic CMG's) — while `ac` keeps a
solve-time edge on the IPM-chimera and several unstructured families. All
three solvers converged to the 1e-8 tolerance on every instance of every
family — not a single failure.

**Worst case per family** (minimum per-instance total-time speedup — how badly
`cmg-k-elim` can lose within each family, and where):

| family | worst vs `ac` | on instance | worst vs `cmg-v` | on instance |
|---|---:|---|---:|---|
| sachdeva_star | 1.02× | star_join(K₁₀₀, 50) | 0.54× | star_join(K₃₀₀, 150) |
| checkered | 2.37× | 2·10⁸ nnz, 2³ blocks | 0.75× | 2·10⁸ nnz, 4³ blocks |
| aniso | 2.03× | 2·10⁸ nnz, ξ=0.001 | 0.92× | 2·10⁸ nnz, ξ=1 |
| uniform_grid | 1.35× | 2·10⁶ nnz | 0.94× | 2·10⁷ nnz |
| wgrid | 0.86× | 2·10⁸ nnz, w=1000 | 0.68× | 2·10⁸ nnz, w=0.01 |
| suitesparse (Laplacian) | 0.83× | Gaertner/nopoly | 0.57× | Andrews/Andrews |
| suitesparse (SDDM) | 0.58× | Oberwolfach/t3dl_e | 0.67× | HB/bcsstm23 |
| chimeraIPM | 0.31× | uc.n1e5.i3.eps1e-5.2 | 0.47× | uc.n1e5.i3.eps1e-5.2 |
| spielmanIPM | 0.84× | k500.i7 | 6.14× | k300.i7 |
| uni_chimera | 0.58× | n=10⁶ | 0.85× | n=10⁶ |
| uni_bndry_chimera | 0.56× | n≈10⁶ | 0.82× | n≈10⁶ |
| wted_chimera | 0.87× | n=10⁶ | 0.85× | n=10⁶ |
| wted_bndry_chimera | 0.69× | n≈10⁴ | 0.83× | n≈10⁷ |

The single worst instance across the whole suite is one chimera-IPM Laplacian
where `cmg-k-elim` takes 3.2× `ac`'s total time (0.31×); on every structured
family its worst instance still ties or beats `ac`, the sole exception being
the heaviest `wgrid` weighting (w=1000, 0.86× — the other six `wgrid`
instances are 1.10–2.40×).

**Solver columns.** `ac` = ApproxChol, the paper's base solver
(`ApproxCholParams(:deg,0,0)`); `ac-s2m2` = its split-2/merge-2 variant (the
paper's AC-s2m2, recommended there for hard problems); `cmg-v` = classic CMG
(stationary cycle + PCG); `cmg-k-elim` = CMG K-cycle with exact degree-1/2
elimination (the default of
[CombinatorialMultigrid.jl](https://github.com/ikoutis/CombinatorialMultigrid.jl)).

### Caveats

- **Environment.** All numbers are from single compute nodes of NJIT's Wulver
  cluster, **single-threaded** (Julia 1.11.9, `JULIA_NUM_THREADS=1`, BLAS pinned
  to 1 thread), seed 1, **3 repetitions**, reporting **medians** — the same
  protocol as the paper (tolerance 1e-8, maxits 1000, `b` random and projected
  to the range space). One machine, one seed: treat small ratios (±10%) as noise.
- **`sachdeva_star` favors CMG for a structural reason** — its clique clusters
  are ideal for CMG's aggregation (constant ~7 iterations at every size). The
  4.35× is measured against the paper's *base* `ac`, whose iteration count blows
  up on this family (it exceeds the 1e-8 target on the largest instances, max
  relres ≈ 7e-6 — the stagnation the paper itself reports for AC here). Against
  the robust `ac-s2m2` the speedup is ≈ 3×.
- **SPE is skipped.** The SPE family needs a manual download we did not stage;
  the coverage report lists it as intentionally skipped, not missing.
- **Chimera rows are per-size medians** over 3 independent random graphs per
  size (the generators draw a fresh graph per repetition).
- This comparison covers the two ApproxChol solvers and the CMG variants only —
  not the paper's full external-solver sweep (LAMG, HyPre, PETSc, ICC, MATLAB-CMG).

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

# 3. Manual data (required for the full IPM sweep shown in the tables):
#    get ipmMat.zip — link in Tutorial.md §IPM — and unzip it in matrix-files/.
#    It extracts into matrix-files/ipmMat/; the loaders read it there directly.
#      cd ../matrix-files && unzip /path/to/ipmMat.zip && cd ../performance-experiments
#    (Optional: SPE needs spe.zip, Tutorial.md §SPE; without it the spe family
#     is skipped and reported as such.)

# 4. Run the benchmark (sequential; pins the four table solvers, 3 reps)
./run_paper_comparison.sh run --scale paper           # or: medium | smoke

# 5. Aggregate into the public tables
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
