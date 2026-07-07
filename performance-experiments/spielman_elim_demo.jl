#==========================================================
Spielman IPM elimination demo.

Loads the real Spielman IPM Laplacians (the `spielmanIPM` family, `sk<k>i<i>`,
downloaded/cached from SuiteSparse) and, for each one, reports the surviving
"core" size after exact degree-1/2 elimination, then compares four solvers on
build time, solve time, total time, iterations, and true relative residual:

    ac          approxchol_lap  (ApproxChol, default)
    ac-s2m2     approxchol_lap  (split=2, merge=2)
    cmg-k       CombinatorialMultigrid K-cycle
    cmg-k-elim  CombinatorialMultigrid K-cycle on the degree-1/2-eliminated core

This reuses the harness's solver factories (cmgSolvers.jl), matrix loaders
(benchFamilies.jl), and per-solver timing (compareSolversCore.jl), so the
numbers are directly comparable to chol_vs_kcycle.jl — with the extra "core%"
column that quantifies how tree-like each matrix is.

Usage (from performance-experiments/, with the project active):
    julia --project=.. spielman_elim_demo.jl
    julia --project=.. spielman_elim_demo.jl --scale medium
    julia --project=.. spielman_elim_demo.jl --names sk100i1,sk200i3,sk300i5
    julia --project=.. spielman_elim_demo.jl --scale paper --limit 10 --maxits 2000

Options (defaults in brackets):
    --scale smoke|medium|paper   which sk<k>i<i> sweep to load        [medium]
    --names a,b,...              explicit matrix names (overrides --scale)
    --limit N                    cap number of matrices
    --solvers a,b,...            subset/superset of the columns
                                 [ac,ac-s2m2,cmg-k,cmg-k-elim]
                                 (also available: cmg-v, cmg-v-elim)
    --tol T                      relative residual tolerance          [1e-8]
    --maxits K                   max (outer) iterations               [1000]
    --seed S                     base seed for RHS + AC randomness    [1]
    --offline                    never download; use only cached data

Prerequisite: the matching CombinatorialMultigrid.jl (with the `eliminate`
branch) must be the active CMG — point setup.jl at it via CMG_REV=main or
CMG_DEV=/path, then run setup.jl before this script.
===========================================================#

using Random
using SparseArrays
using LinearAlgebra
using Statistics
using Printf
using Laplacians
using CombinatorialMultigrid

include(joinpath(@__DIR__, "..", "julia-files", "compareSolversCore.jl"))
include(joinpath(@__DIR__, "..", "julia-files", "cmgSolvers.jl"))
include(joinpath(@__DIR__, "..", "julia-files", "benchFamilies.jl"))

# ----------------------------------------------------------------- arg parsing
function parse_opts(args)
    opts = Dict{Symbol,Any}(
        :scale => :medium, :names => nothing, :limit => nothing,
        :solvers => ["ac", "ac-s2m2", "cmg-k", "cmg-k-elim"],
        :tol => 1e-8, :maxits => 1000, :seed => 1, :offline => false,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        needsval() = (i + 1 <= length(args) ? args[i+1] : error("missing value after $a"))
        if a == "--scale"
            opts[:scale] = Symbol(needsval()); i += 2
        elseif a == "--names"
            opts[:names] = String.(split(needsval(), ',')); i += 2
        elseif a == "--limit"
            opts[:limit] = parse(Int, needsval()); i += 2
        elseif a == "--solvers"
            opts[:solvers] = String.(split(needsval(), ',')); i += 2
        elseif a == "--tol"
            opts[:tol] = parse(Float64, needsval()); i += 2
        elseif a == "--maxits"
            opts[:maxits] = parse(Int, needsval()); i += 2
        elseif a == "--seed"
            opts[:seed] = parse(Int, needsval()); i += 2
        elseif a == "--offline"
            opts[:offline] = true; i += 1
        elseif a == "-h" || a == "--help"
            println("usage: julia --project=.. spielman_elim_demo.jl " *
                    "[--scale smoke|medium|paper] [--names a,b,...] [--limit N] " *
                    "[--solvers ac,ac-s2m2,cmg-k,cmg-k-elim,cmg-v,cmg-v-elim] " *
                    "[--tol T] [--maxits K] [--seed S] [--offline]")
            exit(0)
        else
            error("unknown argument: $a")
        end
    end
    return opts
end

# ------------------------------------------------------------------- instances
function spielman_instances(opts)
    if opts[:names] !== nothing
        insts = [BenchInstance("spielmanIPM $(nm)", (bs, rep) -> begin
                     L = loadIPM(nm)
                     L === nothing && return nothing
                     a, _ = adj(L)
                     (:lap, a, nm)
                 end) for nm in opts[:names]]
        return applyLimit(insts, opts[:limit])
    end
    return spielmanIPMInstances(opts[:scale]; limit = opts[:limit])
end

# ----------------------------------------------------------------------- setup
opts = parse_opts(ARGS)
allow_download!(!opts[:offline])

tests_lap, _ = cholVsKcycleTests(solvers = opts[:solvers])   # matched column set
instances = spielman_instances(opts)

println("solvers: ", join([t.name for t in tests_lap], ", "))
println("scale=$(opts[:scale]) tol=$(opts[:tol]) maxits=$(opts[:maxits]) seed=$(opts[:seed]) offline=$(opts[:offline])")
println("BLAS threads: ", BLAS.get_num_threads(), "\n")

if isempty(instances)
    println("No Spielman matrices available. Fetch them first, e.g.:")
    println("  julia --project=.. performance-experiments/download_data.jl --scale $(opts[:scale])")
    println("or drop sk<k>i<i>.mat/.mm into matrix-files/ and rerun (add --offline).")
    exit(0)
end

# warm up compilation on a tiny near-tree so the first matrix's timings are clean
print("warming up (JIT)… "); flush(stdout)
let
    Random.seed!(0)
    n = 2000
    I = Int[]; J = Int[]; V = Float64[]
    for v = 2:n
        p = rand(1:v-1); w = 0.5 + rand()
        push!(I, v); push!(J, p); push!(V, w); push!(I, p); push!(J, v); push!(V, w)
    end
    aw = sparse(I, J, V, n, n)
    lw = lap(aw); bw = randn(n); bw .-= mean(bw)
    for t in tests_lap
        try
            f = t.solver(aw; tol = 1e-6, maxits = 30, verbose = false)
            f(bw; pcgIts = [0], tol = 1e-6, maxits = 30, verbose = false)
        catch
        end
    end
end
println("done\n")

# --------------------------------------------------------------------- run all
summary = Dict{String,Vector{Float64}}()   # name -> total times, for a closing recap
core_frac = Float64[]

for inst in instances
    loaded = inst.load(opts[:seed], 1)
    loaded === nothing && continue
    kind, a, name = loaded
    la = lap(a)
    n = size(la, 1)
    ne = nnz(la)

    # surviving core after exact degree-1/2 elimination (headline "how tree-like")
    _, ind, _, _ = CombinatorialMultigrid.eliminate_deg12(la)
    core = length(ind)
    push!(core_frac, 100 * core / n)

    Random.seed!(hash((opts[:seed], name)))
    b = randn(n); b .-= mean(b)

    @printf("=== %-10s  n=%-9d  nnz=%-10d  core=%-7d (%.3f%% of n) ===\n",
            name, n, ne, core, 100 * core / n)
    @printf("  %-12s %10s %10s %10s %8s %10s\n",
            "solver", "build_s", "solve_s", "tot_s", "its", "relres")
    for t in tests_lap
        seed = solverSeed(opts[:seed], name, 1, t.name)
        ret = testSolverCore(t.solver, a, la, b, opts[:tol], opts[:maxits], false; seed = seed)
        # ret = (solve_time, build_time, iters, relerr, x)
        tot = ret[1] + ret[2]
        @printf("  %-12s %10.3f %10.3f %10.3f %8d %10.1e\n",
                t.name, ret[2], ret[1], tot, Int(isfinite(ret[3]) ? ret[3] : -1), ret[4])
        push!(get!(summary, t.name, Float64[]), tot)
    end
    println()
end

# ------------------------------------------------------------------- recap
println("Median total time (build+solve) across ", length(core_frac), " matrices:")
for t in tests_lap
    ts = get(summary, t.name, Float64[])
    finite = filter(isfinite, ts)
    med = isempty(finite) ? NaN : median(finite)
    @printf("  %-12s  %8.3f s\n", t.name, med)
end
if !isempty(core_frac)
    @printf("Core fraction: min %.3f%%  median %.3f%%  max %.3f%%\n",
            minimum(core_frac), median(core_frac), maximum(core_frac))
end
println("\nColumns: build_s/solve_s/tot_s in seconds; its = outer iterations")
println("(not work-comparable across ac vs cmg — compare tot_s and relres).")
println("core% = surviving core after exact degree-1/2 elimination; small => near-tree.")
