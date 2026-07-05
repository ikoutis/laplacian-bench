#==========================================================
Approximate-Cholesky vs CMG K-cycle benchmark runner.

Compares the Laplacians.jl approxchol variants ("ac", "ac-s2m2") against
CombinatorialMultigrid's K-cycle ("cmg-k") and V-cycle/PCG ("cmg-v") on any
of the benchmark problem families, with reproducible seeding and repetitions.

    julia --project=.. chol_vs_kcycle.jl <family|all> [options]

Options (defaults in brackets):
    --scale smoke|medium|paper   instance sweep size            [smoke]
    --n 1e6[,1e7,...]            size list for chimera families [per scale]
    --reps R                     repetitions per instance       [3]
    --seed S                     base seed (>= 0)               [1]
    --tol T                      relative residual tolerance    [family default]
    --maxits K                   max (outer) iterations         [1000]
    --solvers a,b,...            subset of ac,ac-s2m2,cmg-k,cmg-v [all]
    --limit N                    cap number of instances per family
    --max-hours H                stop before starting work past H hours
    --out DIR                    output directory [../performance-analyses/chol-vs-kcycle]
    --no-warmup                  skip the warm-up solves

Families: uniform_grid aniso wgrid checkered sachdeva_star suitesparse spe
          chimeraIPM spielmanIPM uni_chimera uni_bndry_chimera wted_chimera
          wted_bndry_chimera

Reproducibility: with the same (family, scale/--n, --reps, --seed, --solvers)
and Julia version, reruns are bit-identical — chimera instances are indexed by
i = 1000*seed + rep, the RHS RNG is seeded per (seed, label, rep), and each
solver's build re-seeds the global RNG per (seed, label, rep, solvername).

Results go to <out>/<name>.jld2 in the same Dict/@save format the existing
performance-analyses notebooks read, saved after every repetition.

Exit codes: 0 success (including nothing-to-do), 1 bad arguments,
2 unexpected runtime failure.
===========================================================#

using Random
using LinearAlgebra
using SparseArrays
using Statistics
using Printf
using Dates
using JLD2
using Laplacians

include(joinpath(@__DIR__, "..", "julia-files", "compareSolversCore.jl"))
include(joinpath(@__DIR__, "..", "julia-files", "cmgSolvers.jl"))
include(joinpath(@__DIR__, "..", "julia-files", "benchFamilies.jl"))

const USAGE = """
usage: julia --project=.. chol_vs_kcycle.jl <family|all>
           [--scale smoke|medium|paper] [--n 1e6[,1e7]] [--reps R] [--seed S]
           [--tol T] [--maxits K] [--solvers ac,ac-s2m2,cmg-k,cmg-v]
           [--limit N] [--max-hours H] [--out DIR] [--no-warmup]
families: $(join(FAMILY_ORDER, " ")), or: all
"""

function parseArgs(args)
    opts = Dict{Symbol,Any}(
        :family => nothing, :scale => :smoke, :n => nothing, :ntoken => nothing,
        :reps => 3, :seed => 1, :tol => nothing, :maxits => 1000,
        :solvers => nothing, :limit => nothing, :maxhours => nothing,
        :out => normpath(joinpath(@__DIR__, "..", "performance-analyses", "chol-vs-kcycle")),
        :warmup => true,
    )
    i = 1
    while i <= length(args)
        a = args[i]
        needsval(flag) = (i += 1) <= length(args) ? args[i] :
                         error("missing value for $(flag)")
        if a == "--scale"
            s = Symbol(needsval(a))
            s in (:smoke, :medium, :paper) || error("bad --scale $(s)")
            opts[:scale] = s
        elseif a == "--n"
            tok = needsval(a)
            opts[:ntoken] = tok
            opts[:n] = [Int64(Float64(Base.parse(Float64, t))) for t in split(tok, ',')]
        elseif a == "--reps"
            opts[:reps] = Base.parse(Int, needsval(a))
        elseif a == "--seed"
            opts[:seed] = Base.parse(Int, needsval(a))
            opts[:seed] >= 0 || error("--seed must be >= 0")
        elseif a == "--tol"
            opts[:tol] = Base.parse(Float64, needsval(a))
        elseif a == "--maxits"
            opts[:maxits] = Base.parse(Int, needsval(a))
        elseif a == "--solvers"
            opts[:solvers] = String.(split(needsval(a), ','))
        elseif a == "--limit"
            opts[:limit] = Base.parse(Int, needsval(a))
        elseif a == "--max-hours"
            opts[:maxhours] = Base.parse(Float64, needsval(a))
        elseif a == "--out"
            opts[:out] = abspath(needsval(a))
        elseif a == "--no-warmup"
            opts[:warmup] = false
        elseif a == "--help" || a == "-h"
            println(USAGE)
            exit(0)
        elseif startswith(a, "--")
            error("unknown option $(a)")
        elseif opts[:family] === nothing
            opts[:family] = a
        else
            error("unexpected argument $(a)")
        end
        i += 1
    end
    opts[:family] === nothing && error("missing <family>")
    opts[:family] == "all" || haskey(FAMILIES, opts[:family]) ||
        error("unknown family $(opts[:family])")
    return opts
end

# One tiny Laplacian + one tiny SDDM instance through every solver, so that
# JIT compilation never lands inside a timed region.
function warmup(tests_lap, tests_sddm)
    println("----- warm up starting ------")
    dicWarmup = Dict()

    a = uni_chimera(1000, 0)
    b = randn(size(a, 1))
    b = lap(a) * b
    b = b ./ norm(b)
    coreTestLap(tests_lap, dicWarmup, a, b;
        tol = 1e-8, maxits = 1000, testName = "warmup_uni_chimera(1000,0)", baseseed = 0)

    M = uniform_grid_sddm(1000)
    bs = randn(size(M, 1))
    bs = M * bs
    bs = bs ./ norm(bs)
    coreTestSddm(tests_sddm, dicWarmup, M, bs;
        tol = 1e-8, maxits = 1000, testName = "warmup_uniform_grid_sddm(1000)", baseseed = 0)

    println("----- warm up complete ------")
end

outName(famname, opts) = begin
    scaletag = opts[:ntoken] !== nothing ? "n" * replace(opts[:ntoken], ',' => '-') :
               String(opts[:scale])
    "$(famname).$(scaletag).seed$(opts[:seed]).reps$(opts[:reps]).jld2"
end

function runFamily(fam::BenchFamily, tests_lap, tests_sddm, opts, t0)
    tol = opts[:tol] === nothing ? fam.tol : opts[:tol]
    seed = opts[:seed]
    reps = opts[:reps]

    if opts[:n] !== nothing && !fam.sized
        @warn "--n ignored for family $(fam.name) (fixed instance sweep)"
    end

    insts = fam.sized ?
        fam.instances(opts[:scale]; n = opts[:n], limit = opts[:limit]) :
        fam.instances(opts[:scale]; limit = opts[:limit])

    if isempty(insts)
        println("family $(fam.name): no instances available at scale $(opts[:scale]) — nothing to do")
        return 0
    end

    fn = joinpath(opts[:out], outName(fam.name, opts))
    println("family $(fam.name): $(length(insts)) instance(s) x $(reps) rep(s), tol=$(tol)")
    println("results -> $(fn)")

    dic = Dict()
    dic["run_args"] = join(ARGS, " ")
    dic["julia_version"] = string(VERSION)
    ran = 0

    for inst in insts
        for rep in 1:reps
            if opts[:maxhours] !== nothing && (time() - t0) > 3600 * opts[:maxhours]
                println("--max-hours budget reached; stopping with partial results")
                return ran
            end

            # Generators are deterministic today, but pin the RNG anyway in
            # case a future Laplacians version draws from it during generation.
            Random.seed!(solverSeed(seed, inst.base, rep, "instance"))
            loaded = inst.load(seed, rep)
            if loaded === nothing
                @warn "skipping unavailable instance" family = fam.name instance = inst.base
                break   # data won't appear between reps
            end
            kind, mat, label = loaded
            sys = kind === :lap ? lap(mat) : mat

            println("=============================================")
            println("$(fam.name) | $(label) | rep $(rep)/$(reps) | nv=$(size(sys,1)) nnz=$(nnz(sys))")

            Random.seed!(solverSeed(seed, label, rep, "rhs"))
            b = sys * randn(size(sys, 1))
            b = b ./ norm(b)

            if kind === :lap
                coreTestLap(tests_lap, dic, mat, b;
                    tol = tol, maxits = opts[:maxits], testName = label,
                    rep = rep, baseseed = seed, la = sys)
            else
                coreTestSddm(tests_sddm, dic, mat, b;
                    tol = tol, maxits = opts[:maxits], testName = label,
                    rep = rep, baseseed = seed)
            end

            dic["saved_at"] = string(Dates.now())
            @save fn dic
            ran += 1

            # Large instances: drop references before generating the next one.
            mat = nothing; sys = nothing; b = nothing
            GC.gc()
        end
    end

    return ran
end

function main()
    opts = try
        parseArgs(ARGS)
    catch e
        println(stderr, "argument error: ", sprint(showerror, e))
        println(stderr, USAGE)
        exit(1)
    end

    tests_lap, tests_sddm = try
        cholVsKcycleTests(solvers = opts[:solvers])
    catch e
        println(stderr, "argument error: ", sprint(showerror, e))
        exit(1)
    end
    println("solvers: ", join([t.name for t in tests_lap], ", "))
    println("scale=$(opts[:scale]) reps=$(opts[:reps]) seed=$(opts[:seed]) maxits=$(opts[:maxits])")

    mkpath(opts[:out])

    t0 = time()
    opts[:warmup] && warmup(tests_lap, tests_sddm)

    famnames = opts[:family] == "all" ? FAMILY_ORDER : [opts[:family]]
    total = 0
    for famname in famnames
        total += runFamily(FAMILIES[famname], tests_lap, tests_sddm, opts, t0)
    end

    @printf("done: %d benchmark run(s) in %.1f minutes\n", total, (time() - t0) / 60)
end

try
    main()
catch e
    e isa InterruptException && rethrow()
    println(stderr, "FATAL: ", sprint(showerror, e, catch_backtrace()))
    exit(2)
end
