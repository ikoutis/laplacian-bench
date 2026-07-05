#==========================================================
Slim, pure-Julia benchmark harness.

This is the native-solver core of compareSolvers.jl, without the
MATLAB/PETSc/HyPre bridges (which drag in MATLAB.jl, CSV, DataFrames and
external binaries). It exists so the chol-vs-kcycle experiments can run on
machines with nothing but Julia — e.g. Wulver compute nodes.

Differences from compareSolvers.jl, all additive:
  * no external-solver keywords (test_cmg, test_icc, ...);
  * optional deterministic seeding: the global RNG is re-seeded right before
    each solver's build, derived from (baseseed, testName, rep, solvername),
    so randomized solvers (approxchol) are reproducible per column and results
    for one column don't shift when another column is added or removed.
    Reproducibility is per Julia version (Base.hash is not stable across
    versions);
  * two extra bookkeeping columns, "rep" and "seed" (old notebooks index
    specific keys, so they are unaffected);
  * function names differ (testSolverCore, coreTestLap, coreTestSddm) so this
    file can be loaded alongside compareSolvers.jl without method clashes.

The result-dictionary schema is otherwise identical:
  nv, ne, testName, names, and {name}_solve/_build/_tot/_its/_err.
===========================================================#

using Laplacians
using SparseArrays
using LinearAlgebra
using Statistics
using Random

"""
    SolverTest(solver, name)

Encloses a solver with its name, so that we can compare it in tests.
The solver contract (same as compareSolvers.jl / Laplacians.jl):
`f = solver(a; tol, maxits, verbose)` builds, then `x = f(b; pcgIts, tol,
maxits, verbose)` solves, writing the iteration count into `pcgIts[1]`.
"""
struct SolverTest
    solver::Function
    name::String
end

"""
    initDictCol!(dic, name, typ)

For a dictionary in which each key indexes an array.
If dic does not contain an entry of `name`, create with set to `Array(typ,0)`.
"""
function initDictCol!(dic, name, typ)
    if ~haskey(dic, name)
        dic[name] = typ[]
    end
end

"""
`ret` is the answer returned by a speed test.
This pushes it into the dictionary on which we are storing the tests.
"""
function pushSpeedResult!(dic, name, ret)
    push!(dic["$(name)_solve"], ret[1])
    push!(dic["$(name)_build"], ret[2])
    push!(dic["$(name)_tot"], ret[1] + ret[2])
    push!(dic["$(name)_its"], ret[3])
    push!(dic["$(name)_err"], ret[4])
end

# Deterministic per-solver seed. `baseseed === nothing` disables seeding.
solverSeed(baseseed, testName, rep, name) =
    hash((UInt64(baseseed), String(testName), Int(rep), String(name)))

#=
Time one native solver on one system. `mat` is whatever the solver consumes
(an adjacency matrix for Laplacian solvers, the SDDM matrix itself for SDDM
solvers); `sys` is the system matrix used for the residual check (lap(mat) in
the Laplacian case, mat itself in the SDDM case). Returns
(solve_time, build_time, iters, relerr, x), or (Inf,Inf,Inf,Inf,Inf) if the
solver threw. This is a general safety net: the curated benchmark matrices are
all SDDM / near-SDDM (see Tutorial.md), so no benchmark input trips CMG's
validateInput!, but it guards against an asymmetric or arbitrary user-supplied
matrix (validateInput! throws on asymmetry or a positive off-diagonal).
=#
function testSolverCore(solver, mat, sys, b, tol, maxits, verbose; seed = nothing)

    try
        GC.gc()
        if seed !== nothing
            Random.seed!(seed)
        end
        t0 = time()
        f = solver(mat, tol = tol, maxits = maxits, verbose = verbose)
        build_time = time() - t0

        it = [0]
        GC.gc()

        t0 = time()
        x = f(b, pcgIts = it, tol = tol, maxits = maxits, verbose = verbose)
        solve_time = time() - t0

        err = norm(sys * x .- b) / norm(b)

        ret = (solve_time, build_time, it[1], err, x)
        if verbose
            println("Solve time, build time, iter, err:", (solve_time, build_time, it[1], err))
        end
        return ret
    catch e
        println("Solver Error: ", sprint(showerror, e))
        return (Inf, Inf, Inf, Inf, Inf)
    end

end

# Shared body of coreTestLap / coreTestSddm.
function coreTestRun(solvers, dic::Dict, mat, sys, b;
        tol = 1e-8, maxits = 1000, verbose = false, testName = "",
        rep = 1, baseseed = nothing)

    initDictCol!(dic, "nv", Int)
    initDictCol!(dic, "ne", Int)
    initDictCol!(dic, "testName", String)
    initDictCol!(dic, "rep", Int)
    initDictCol!(dic, "seed", Int)

    solvecol(name) = "$(name)_solve"
    buildcol(name) = "$(name)_build"
    totcol(name) = "$(name)_tot"
    itscol(name) = "$(name)_its"
    errcol(name) = "$(name)_err"

    dic["names"] = String[]
    for t in solvers
        push!(dic["names"], t.name)
    end

    for name in dic["names"]
        initDictCol!(dic, solvecol(name), Float64)
        initDictCol!(dic, buildcol(name), Float64)
        initDictCol!(dic, totcol(name), Float64)
        initDictCol!(dic, itscol(name), Float64)
        initDictCol!(dic, errcol(name), Float64)
    end

    push!(dic["nv"], size(sys, 1))
    push!(dic["ne"], nnz(sys))
    push!(dic["testName"], testName)
    push!(dic["rep"], Int(rep))
    push!(dic["seed"], baseseed === nothing ? -1 : Int(baseseed))

    x = []

    for i in 1:length(solvers)
        solverTest = solvers[i]

        println("--------------")
        println(solverTest.name)

        seed = baseseed === nothing ? nothing :
               solverSeed(baseseed, testName, rep, solverTest.name)
        ret = testSolverCore(solverTest.solver, mat, sys, b, tol, maxits, verbose;
                             seed = seed)

        if i == 1
            x = ret[5]
        end
        println("total: $(ret[1]+ret[2]), iter: $(ret[3]), solve: $(ret[1]), build: $(ret[2]), err: $(ret[4])")
        println("--------------")
        pushSpeedResult!(dic, solverTest.name, ret)
    end

    return x
end

"""
    coreTestLap(solvers, dic, a, b; tol, maxits, verbose, testName, rep, baseseed, la)

Native-solver-only counterpart of `testLap`: `a` is an adjacency matrix, each
solver receives `a`, and residuals are computed against `lap(a)`. `b` is
mean-centered (Laplacian systems are only solvable on the range space).
Pass `la` to reuse an already-computed `lap(a)` (worthwhile at 1e7 nodes).
Results accumulate into `dic` (same schema as compareSolvers.jl); returns the
first solver's solution.
"""
function coreTestLap(solvers, dic::Dict, a::SparseMatrixCSC, b::Array;
        tol = 1e-8, maxits = 1000, verbose = false, testName = "",
        rep = 1, baseseed = nothing, la = nothing)

    b = b .- mean(b)
    if la === nothing
        la = Laplacians.lap(a)
    end

    coreTestRun(solvers, dic, a, la, b;
        tol = tol, maxits = maxits, verbose = verbose, testName = testName,
        rep = rep, baseseed = baseseed)
end

"""
    coreTestSddm(solvers, dic, sddmmat, b; tol, maxits, verbose, testName, rep, baseseed)

Native-solver-only counterpart of `testSddm`: each solver receives the SDDM
matrix itself, residuals are computed against it, and `b` is used as given.
"""
function coreTestSddm(solvers, dic::Dict, sddmmat::SparseMatrixCSC, b::Array;
        tol = 1e-8, maxits = 1000, verbose = false, testName = "",
        rep = 1, baseseed = nothing)

    coreTestRun(solvers, dic, sddmmat, sddmmat, b;
        tol = tol, maxits = maxits, verbose = verbose, testName = testName,
        rep = rep, baseseed = baseseed)
end
