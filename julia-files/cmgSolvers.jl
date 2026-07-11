#==========================================================
SolverTest factories for the chol-vs-kcycle comparison:

  * approxchol (Laplacians.jl)            -> columns "ac", "ac-s2m2", ...
  * CombinatorialMultigrid K-cycle        -> column  "cmg-k"
  * CombinatorialMultigrid V-cycle (PCG)  -> column  "cmg-v"

Both CMG columns are driven by `cmg_solve`, whose flexible-CG outer loop
reduces to plain PCG for `:vcycle`, so the two cycles differ only in the
preconditioning cycle — the fair comparison of what the K-cycle buys. The
K-cycle preconditioner is NONLINEAR, so it must never be used inside a
standard PCG (e.g. Laplacians.pcgSolver); `cmg_solve` is the only driver.

The hierarchy is built once inside the harness's timed build phase (via the
public `cmg_preconditioner_lap`, which validates the input and, for strictly
dominant SDD matrices, grounds them with an extra coordinate) and reused by
`cmg_solve(H, b)`, which handles the padding/extraction itself — so `b` is
always passed at its original size.

Column names deliberately avoid "cmg"/"cmg2", which older notebooks use for
the MATLAB CMG solver.

Iteration-count caveat: `cmg-k_its` counts OUTER flexible-CG iterations, each
of which hides budget-capped inner FCG work at the coarse levels. Iteration
counts are not work-comparable across columns; compare `_tot`/`_solve`, and
use `_err` (computed uniformly by the harness) as the accuracy metric.
===========================================================#

using Laplacians
using CombinatorialMultigrid
using SparseArrays
using LinearAlgebra

cmg_cycle_name(cycle::Symbol) = cycle === :kcycle ? "cmg-k" :
                                cycle === :vcycle ? "cmg-v" :
                                error("unknown CMG cycle $(repr(cycle))")

# ---------------------------------------------------------------------------
# Injection counting for the sparsify-on-stall columns.
#
# A "sparsify" level is a same-size hierarchy level (identity transfer,
# `cI = 1:n`, `nc == n`) injected when aggregation stalls. It is detected
# exactly as the package does: `!islast && nc == n`. `n_inject` walks whatever
# hierarchy object `cmg_preconditioner_lap` returned — a bare
# `Vector{HierarchyLevel}` (:kcycle/:vcycle path), an `EliminatedHierarchy`
# (`eliminate=true`), or a `DisconnectedHierarchy` (`split_components`, whose
# blocks are per-component hierarchies or a scalar for a singleton) — and sums
# the injected levels. Duck-typed on fields so it needs no non-exported types.
function n_inject(H)::Int
    if H isa AbstractVector                       # Vector{HierarchyLevel}
        return count(h -> (!h.islast && h.nc == h.n), H)
    elseif hasproperty(H, :H)                     # EliminatedHierarchy
        return n_inject(H.H)
    elseif hasproperty(H, :blocks)                # DisconnectedHierarchy
        return sum(n_inject(b) for b in H.blocks; init = 0)
    else                                          # scalar singleton block, etc.
        return 0
    end
end

# Callable wrapper around a CMG solve closure that also advertises the number
# of sparsify levels injected at build time. The harness (testSolverCore) reads
# `f.n_inject` if present, so `ac`/other solvers (plain closures) record 0.
struct CMGSolve{F}
    f::F
    n_inject::Int
end
(s::CMGSolve)(b; kwargs...) = s.f(b; kwargs...)

# "-elim" columns first exactly factor out degree-1/2 nodes (partial Cholesky)
# and run the chosen cycle on the surviving core; see cmg_preconditioner_lap's
# `eliminate` keyword. Especially effective on near-tree graphs (e.g. the
# Spielman IPM family) where the core is tiny.
cmg_elim_name(cycle::Symbol) = cmg_cycle_name(cycle) * "-elim"

# Shared solve-closure constructor: `sys` is the SDD system CMG works on
# (lap(a) for the Laplacian path, M itself for the SDDM path). With
# `eliminate = true`, `cmg_preconditioner_lap` returns an EliminatedHierarchy
# and `cmg_solve` dispatches on it — `b` is still passed at its original size.
function cmg_build(sys, cycle, theta, inner_tol, tol, maxits; eliminate = false,
                  sparsify_on_stall = false, sparsify_opts = nothing)
    theta, inner_tol = Float64(theta), Float64(inner_tol)
    # Non-sparsify columns call cmg_preconditioner_lap exactly as before (no
    # sparsify kwargs, no SparsifyOptions reference) so they are unaffected by
    # this change. Only the sparsify columns pass the sparsify options.
    (_, H) = if sparsify_on_stall
        cmg_preconditioner_lap(sys; cycle = cycle, theta = theta,
            inner_tol = inner_tol, eliminate = eliminate, sparsify_on_stall = true,
            sparsify_opts = sparsify_opts === nothing ? SparsifyOptions() : sparsify_opts)
    else
        cmg_preconditioner_lap(sys; cycle = cycle, theta = theta,
            inner_tol = inner_tol, eliminate = eliminate)
    end
    solve = function (b; pcgIts = Int[0], tol = tol, maxits = maxits, verbose = false, args...)
        x, stats = cmg_solve(H, b; tol = Float64(tol), maxit = Int64(maxits),
                             cycle = cycle, theta = theta, inner_tol = inner_tol)
        if length(pcgIts) > 0
            pcgIts[1] = stats.iterations
        end
        if verbose || !stats.converged
            println("CMG $(cycle)$(eliminate ? "-elim" : "")$(sparsify_on_stall ? "+sparsify" : ""): its=$(stats.iterations) relres=$(stats.relres) converged=$(stats.converged)")
        end
        return x
    end
    return CMGSolve(solve, n_inject(H))
end

"""
    make_cmg_lap(cycle; theta=0.75, inner_tol=0.25)

Solver factory for the Laplacian path: receives an ADJACENCY matrix `a`
(the `testLap`/`coreTestLap` convention), builds the CMG hierarchy on
`lap(a)`, and solves with `cmg_solve`.
"""
function make_cmg_lap(cycle::Symbol; theta = 0.75, inner_tol = 0.25, eliminate = false)
    return function (a; tol = 1e-8, maxits = 1000, verbose = false, args...)
        cmg_build(lap(a), cycle, theta, inner_tol, tol, maxits; eliminate = eliminate)
    end
end

"""
    make_cmg_sddm(cycle; theta=0.75, inner_tol=0.25)

Solver factory for the SDDM path: receives the SDDM matrix itself.
`cmg_preconditioner_lap` accepts SDD matrices directly (strictly dominant
rows are grounded via an augmented coordinate). The benchmark's SuiteSparse
selection is all SDDM / near-SDDM with non-positive off-diagonals (see
Tutorial.md), so CMG accepts every benchmark input; validateInput! only
rejects an asymmetric matrix or a positive off-diagonal, and the harness
try/catch would record such a case as an (Inf,...) row.
"""
function make_cmg_sddm(cycle::Symbol; theta = 0.75, inner_tol = 0.25, eliminate = false)
    return function (M; tol = 1e-8, maxits = 1000, verbose = false, args...)
        cmg_build(M, cycle, theta, inner_tol, tol, maxits; eliminate = eliminate)
    end
end

# ---------------------------------------------------------------------------
# Sparsify-on-stall columns. When aggregation stalls (a level that will not
# coarsen because contraction fill densified it, or an expander-like core), the
# build injects a spanner + uniform spectral-sparsifier as a same-size level so
# standard aggregation can resume. Two solve cycles are compared:
#   :legacy  -> "cmg-sparsify-l"  (stationary + PCG, the empirically best cycle)
#   :kscycle -> "cmg-sparsify-ks" (sparsify-aware K-cycle: extra inner FCG at
#                                   injected same-size levels)
# Both build with `eliminate=true` (sparsify + degree-1/2 elimination) per the
# benchmark design. On instances that never stall these reduce to the ordinary
# elimination hierarchy (0 injections) and match the `cmg-k-elim` reference.
cmg_sparsify_name(cycle::Symbol) = cycle === :legacy  ? "cmg-sparsify-l"  :
                                   cycle === :kscycle ? "cmg-sparsify-ks" :
                                   error("unknown sparsify cycle $(repr(cycle))")

# `sparsify_opts = nothing` -> cmg_build constructs the default SparsifyOptions()
# (only on the sparsify path); pass an explicit SparsifyOptions to tune knobs.
function make_cmg_sparsify_lap(cycle::Symbol; theta = 0.75, inner_tol = 0.25,
                               eliminate = true, sparsify_opts = nothing)
    return function (a; tol = 1e-8, maxits = 1000, verbose = false, args...)
        cmg_build(lap(a), cycle, theta, inner_tol, tol, maxits;
                  eliminate = eliminate, sparsify_on_stall = true,
                  sparsify_opts = sparsify_opts)
    end
end

function make_cmg_sparsify_sddm(cycle::Symbol; theta = 0.75, inner_tol = 0.25,
                                eliminate = true, sparsify_opts = nothing)
    return function (M; tol = 1e-8, maxits = 1000, verbose = false, args...)
        cmg_build(M, cycle, theta, inner_tol, tol, maxits;
                  eliminate = eliminate, sparsify_on_stall = true,
                  sparsify_opts = sparsify_opts)
    end
end

# approxchol naming rule, exactly as in the *_ac.jl scripts.
ac_name(split, merge) = split >= 1 && merge >= 1 ? "ac-s$(split)m$(merge)" :
                        split >= 1               ? "ac-s$(split)"          : "ac"

function make_ac_lap(split, merge)
    return function (a; verbose = false, args...)
        approxchol_lap(a; params = ApproxCholParams(:deg, split, merge),
                       verbose = verbose, args...)
    end
end

function make_ac_sddm(split, merge)
    return function (M; verbose = false, args...)
        approxchol_sddm(M; params = ApproxCholParams(:deg, split, merge),
                        verbose = verbose, args...)
    end
end

"""
    (tests_lap, tests_sddm) = cholVsKcycleTests(; ac_pairs=[(0,0),(2,2)],
                                                  cycles=[:kcycle,:vcycle],
                                                  theta=0.75, inner_tol=0.25,
                                                  solvers=nothing)

Build matched Laplacian/SDDM SolverTest arrays: the approxchol variants first
(the first solver's solution is the harness's reference `x`), then the CMG
columns, then the degree-1/2-elimination CMG columns (`cmg-k-elim`,
`cmg-v-elim`), then the sparsify-on-stall columns (`cmg-sparsify-l`,
`cmg-sparsify-ks`). `solvers`, if given, is a list of column names to keep, e.g.
["ac", "cmg-k-elim", "cmg-sparsify-l", "cmg-sparsify-ks"]. `elim_cycles` selects
which cycles get an elimination column (default: same as `cycles`; pass `[]` to
omit them); `sparsify_cycles` selects the sparsify columns (default
`[:legacy, :kscycle]`; pass `[]` to omit them).
"""
function cholVsKcycleTests(; ac_pairs = [(0, 0), (2, 2)],
                             cycles = [:kcycle, :vcycle],
                             elim_cycles = [:kcycle, :vcycle],
                             sparsify_cycles = [:legacy, :kscycle],
                             theta = 0.75, inner_tol = 0.25,
                             solvers = nothing)
    tests_lap = SolverTest[]
    tests_sddm = SolverTest[]

    for (s, m) in ac_pairs
        push!(tests_lap, SolverTest(make_ac_lap(s, m), ac_name(s, m)))
        push!(tests_sddm, SolverTest(make_ac_sddm(s, m), ac_name(s, m)))
    end
    for cyc in cycles
        push!(tests_lap, SolverTest(make_cmg_lap(cyc; theta = theta, inner_tol = inner_tol),
                                    cmg_cycle_name(cyc)))
        push!(tests_sddm, SolverTest(make_cmg_sddm(cyc; theta = theta, inner_tol = inner_tol),
                                     cmg_cycle_name(cyc)))
    end
    for cyc in elim_cycles
        push!(tests_lap, SolverTest(make_cmg_lap(cyc; theta = theta, inner_tol = inner_tol,
                                                 eliminate = true), cmg_elim_name(cyc)))
        push!(tests_sddm, SolverTest(make_cmg_sddm(cyc; theta = theta, inner_tol = inner_tol,
                                                   eliminate = true), cmg_elim_name(cyc)))
    end
    for cyc in sparsify_cycles
        push!(tests_lap, SolverTest(make_cmg_sparsify_lap(cyc; theta = theta,
                                    inner_tol = inner_tol), cmg_sparsify_name(cyc)))
        push!(tests_sddm, SolverTest(make_cmg_sparsify_sddm(cyc; theta = theta,
                                     inner_tol = inner_tol), cmg_sparsify_name(cyc)))
    end

    if solvers !== nothing
        keep = Set(String.(solvers))
        found = Set(t.name for t in tests_lap)
        unknown = setdiff(keep, found)
        isempty(unknown) || error("unknown solver name(s): $(join(sort!(collect(unknown)), ", ")); " *
                                  "available: $(join([t.name for t in tests_lap], ", "))")
        filter!(t -> t.name in keep, tests_lap)
        filter!(t -> t.name in keep, tests_sddm)
    end

    return tests_lap, tests_sddm
end
