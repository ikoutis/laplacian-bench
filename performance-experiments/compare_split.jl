#!/usr/bin/env julia
# Measure what CombinatorialMultigrid's `split_components` knob buys on the
# disconnected chimera draws that defeat the connected-graph solver.
#
# For each difficult draw it rebuilds the exact graph (via the benchmark's own
# generators) and times build + solve with split_components = true vs false,
# reporting iterations / relres / converged. With the knob OFF these are the
# cases from the paper run that error at build (isolated vertex) or stall at the
# iteration cap (multiple components); with it ON each connected component is
# solved independently. Requires a CombinatorialMultigrid that has the knob
# (branch claude/cmg-degree-elimination-qwi3iq).
#
# Usage (compute node, after `source wulver/env_wulver.sh`):
#   julia --project=.. compare_split.jl [--maxit K] [--tol T]
#         [--case fam:n:i]...        # override the default difficult set

using CombinatorialMultigrid, SparseArrays, LinearAlgebra, Random, Printf
include(joinpath(@__DIR__, "..", "julia-files", "benchFamilies.jl"))  # generators

# family => (generator, is-Laplacian-adjacency)
const GEN = Dict(
    "uni_chimera"        => (uni_chimera, true),
    "wted_chimera"       => (wted_chimera, true),
    "uni_bndry_chimera"  => (uni_bndry_chimera_fixed, false),
    "wted_bndry_chimera" => (wted_bndry_chimera_fixed, false),
)

# Difficult disconnected draws found in the paper-count run (all n=1e5):
#   uni_chimera i=35        -> isolated vertex, build error when split=false
#   uni_chimera i=34/72     -> 2 components, 1000-iter stall when split=false
#   uni_bndry_chimera i=72  -> 898 components
#   wted_chimera i=66       -> 2 components
const DEFAULT_CASES = [
    ("uni_chimera", 100000, 35), ("uni_chimera", 100000, 34),
    ("uni_chimera", 100000, 72), ("uni_bndry_chimera", 100000, 72),
    ("wted_chimera", 100000, 66),
]

maxit, tol, cases = 1000, 1e-8, Tuple{String,Int,Int}[]
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if     a == "--maxit"; global maxit = parse(Int, ARGS[i+1]); i += 2
        elseif a == "--tol";   global tol = parse(Float64, ARGS[i+1]); i += 2
        elseif a == "--case"
            p = split(ARGS[i+1], ':'); push!(cases, (String(p[1]), parse(Int, p[2]), parse(Int, p[3]))); i += 2
        else; error("unknown arg $(a)")
        end
    end
end
isempty(cases) && (cases = DEFAULT_CASES)

function run_case(fam, n, k)
    local gen, islap = GEN[fam]
    local sys = islap ? lap(gen(n, k)) : gen(n, k)
    local nv = size(sys, 1)
    Random.seed!(1234)
    local b = sys * randn(nv); b ./= norm(b)
    local ncomp = maximum(CombinatorialMultigrid.components(dropzeros(sys)))
    @printf("== %s(%d,%d)  nv=%d  components=%d ==\n", fam, n, k, nv, ncomp)
    for split in (true, false)
        try
            local H
            local tb = @elapsed ((_, H) = cmg_preconditioner_lap(sys; cycle = :kcycle,
                eliminate = true, split_components = split))
            local x, st
            local ts = @elapsed ((x, st) = cmg_solve(H, b; cycle = :kcycle, maxit = maxit, tol = tol))
            @printf("   split=%-5s  build=%.3g  solve=%.3g  tot=%.3g  its=%d  relres=%.2e  converged=%s\n",
                split, tb, ts, tb + ts, st.iterations, st.relres, st.converged)
        catch e
            @printf("   split=%-5s  ERROR: %s\n", split, sprint(showerror, e))
        end
    end
end

# JIT warm-up on a tiny disconnected graph so the first timed case is clean.
let
    W = blockdiag(lap(uni_chimera(1000, 0)), spzeros(1, 1))
    bw = W * randn(size(W, 1))
    for s in (true, false)
        try cmg_solve(W, bw; split_components = s, eliminate = true, maxit = 20) catch end
    end
end

for (f, n, k) in cases
    run_case(f, n, k)
end
