#!/usr/bin/env julia
# Reconstruct and characterize the chimera draws where CMG converged but was much
# slower than ac. Auto-finds the worst converged cmg-k-elim instances per
# (family, size) from the raw .jld2, rebuilds each exact graph via the benchmark's
# own generators, and reports:
#   - structure: components, degree profile, deg<=2 fraction, weight ratio, hubs;
#   - the CMG coarsening hierarchy (level sizes + ratios) — poor coarsening is the
#     usual cause of a weak preconditioner / many iterations;
#   - the recorded ac vs cmg-k / cmg-k-elim build/solve/its breakdown, so the
#     slowness is attributed to build vs iterations.
#
# Usage (compute node, after `source wulver/env_wulver.sh`):
#   julia --project=.. diagnose_chimera_slow.jl [RESULTS_DIR] [--seed S] [--reps R]
#                                               [--top K] [--max-n N] [--only FAM]
# Defaults: RESULTS_DIR=../performance-analyses/chol-vs-kcycle, seed 2, reps 1,
#           top 2 per size, max-n 1_000_000 (skip 1e7 by default), all families.

using JLD2, SparseArrays, LinearAlgebra, Printf
include(joinpath(@__DIR__, "..", "julia-files", "benchFamilies.jl"))  # generators + Laplacians
using CombinatorialMultigrid

const SOLVERS = ["ac", "cmg-k", "cmg-k-elim"]
const FIELDS  = ["build", "solve", "tot", "its"]
# family => (generator, is-laplacian-adjacency)
const GEN = Dict(
    "uni_chimera"        => (uni_chimera, true),
    "wted_chimera"       => (wted_chimera, true),
    "uni_bndry_chimera"  => (uni_bndry_chimera_fixed, false),
    "wted_bndry_chimera" => (wted_bndry_chimera_fixed, false),
)

dir = "../performance-analyses/chol-vs-kcycle"
seed, reps, topk, maxn, only = 2, 1, 2, 1_000_000, ""
let i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if     a == "--seed";  global seed = parse(Int, ARGS[i+1]); i += 2
        elseif a == "--reps";  global reps = parse(Int, ARGS[i+1]); i += 2
        elseif a == "--top";   global topk = parse(Int, ARGS[i+1]); i += 2
        elseif a == "--max-n"; global maxn = parse(Int, ARGS[i+1]); i += 2
        elseif a == "--only";  global only = ARGS[i+1]; i += 2
        else;  global dir = a; i += 1
        end
    end
end

function pool(fam)
    pat = Regex("^" * fam * "\\..*seed$(seed)\\.reps$(reps)\\.jld2\$")
    bysize = Dict{Int,Vector{Any}}()
    for f in sort(readdir(dir; join = true))
        occursin(pat, basename(f)) || continue
        dic = try JLD2.load(f)["dic"] catch; continue end
        haskey(dic, "testName") || continue
        names = get(dic, "names", SOLVERS)
        for (r, tn) in enumerate(dic["testName"])
            nv = haskey(dic, "nv") ? Int(dic["nv"][r]) : -1
            m = Dict{String,Any}()
            for s in SOLVERS
                s in names || continue
                m[s] = Dict(fld => Float64(dic["$(s)_$(fld)"][r]) for fld in FIELDS)
            end
            push!(get!(bysize, nv, Any[]), (name = String(tn), s = m))
        end
    end
    return bysize
end

conv(row, s) = haskey(row.s, s) && isfinite(row.s[s]["tot"]) && row.s[s]["tot"] > 0
ratio(row) = conv(row, "ac") && conv(row, "cmg-k-elim") ?
             row.s["ac"]["tot"] / row.s["cmg-k-elim"]["tot"] : Inf

# adjacency (positive weights) from the system matrix's off-diagonal
function adjacency(sys)
    off = sys - spdiagm(0 => diag(sys))
    a = -off
    dropzeros!(a)
    return a
end

function component_sizes(a)
    n = size(a, 1); seen = falses(n); rows = rowvals(a); sizes = Int[]
    for s in 1:n
        seen[s] && continue
        cnt = 0; stack = [s]; seen[s] = true
        while !isempty(stack)
            u = pop!(stack); cnt += 1
            for idx in nzrange(a, u)
                v = rows[idx]
                if !seen[v]; seen[v] = true; push!(stack, v); end
            end
        end
        push!(sizes, cnt)
    end
    return sort(sizes; rev = true)
end

function structure(a)
    n = size(a, 1); deg = diff(a.colptr); w = nonzeros(a); comps = component_sizes(a)
    @printf("     structure: n=%d m=%d avg_deg=%.2f  deg(min=%d max=%d)  deg<=2=%.1f%%\n",
            n, nnz(a) ÷ 2, sum(deg) / n, minimum(deg), maximum(deg),
            100 * count(<=(2), deg) / n)
    @printf("     weights: min=%.3g max=%.3g ratio=%.3g   components: %d (largest %s)\n",
            minimum(w), maximum(w), maximum(w) / minimum(w),
            length(comps), join(comps[1:min(6, end)], ","))
    p = partialsortperm(deg, 1:min(6, n); rev = true)
    @printf("     top degrees: %s\n", join(deg[p], ","))
end

function coarsening(sys)
    try
        (_, H) = cmg_preconditioner_lap(sys; cycle = :kcycle, eliminate = false)
        ns = [Int(h.n) for h in H]
        rs = [@sprintf("%.2f", ns[i+1] / ns[i]) for i in 1:length(ns)-1]
        @printf("     CMG hierarchy: %d levels  sizes=%s\n", length(ns), join(ns, "->"))
        @printf("     coarsening ratios: %s\n", join(rs, ","))
    catch e
        @printf("     CMG hierarchy: BUILD ERROR: %s\n", sprint(showerror, e))
    end
end

breakdown(row) = for s in SOLVERS
    haskey(row.s, s) || continue
    d = row.s[s]; bs = isfinite(d["tot"]) && d["tot"] > 0 ? 100 * d["build"] / d["tot"] : NaN
    @printf("     %-11s build=%.3g solve=%.3g tot=%.3g its=%.0f (build %.0f%%)\n",
            s, d["build"], d["solve"], d["tot"], d["its"], bs)
end

for fam in sort(collect(keys(GEN)))
    (only == "" || only == fam) || continue
    gen, islap = GEN[fam]
    bysize = pool(fam)
    isempty(bysize) && (println("\n$(fam): no seed$(seed).reps$(reps) files in $(dir)"); continue)
    println("\n=== $(fam) ===")
    for nv in sort(collect(keys(bysize)))
        nv > maxn && (println("  n=$(nv): skipped (> --max-n $(maxn))"); continue)
        cr = sort(filter(r -> conv(r, "cmg-k-elim"), bysize[nv]); by = ratio)
        for row in cr[1:min(topk, end)]
            mm = match(r"\((\d+),\s*(\d+)\)", row.name)
            mm === nothing && (println("  ? cannot parse $(row.name)"); continue)
            n, k = parse(Int, mm.captures[1]), parse(Int, mm.captures[2])
            @printf("  %-26s  ac/cmg-k-elim tot = %.3fx\n", row.name, ratio(row))
            sys = islap ? lap(gen(n, k)) : gen(n, k)
            structure(adjacency(sys))
            coarsening(sys)
            breakdown(row)
        end
    end
end
