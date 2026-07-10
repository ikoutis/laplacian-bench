#!/usr/bin/env julia
# Per-sample tail analysis for the chimera families.
#
# The paper-table CSV collapses each chimera size to one median row and drops
# non-converged samples, so it hides the tail. This reads the raw per-sample
# .jld2 files (including the chunked ones) and, per (family, size), reports:
#   - sample count and per-solver non-convergence count (tot = Inf),
#   - among converged samples, the min / median / max of ac_tot / solver_tot
#     (min = the worst case for the CMG solver),
#   - for each non-converged CMG sample, whether ac converged on that same graph.
#
# Usage (on a compute node, after `source wulver/env_wulver.sh`):
#   julia --project=.. analyze_chimera_tail.jl [RESULTS_DIR] [--seed S] [--reps R]
# Defaults: RESULTS_DIR=../performance-analyses/chol-vs-kcycle, seed 2, reps 1.

using JLD2, Statistics, Printf

const CHIMERA = ["uni_chimera", "uni_bndry_chimera", "wted_chimera", "wted_bndry_chimera"]
const SOLVERS = ["ac", "ac-s2m2", "cmg-v", "cmg-k", "cmg-k-elim"]

dir = "../performance-analyses/chol-vs-kcycle"
seed, reps = 2, 1
i = 1
args = ARGS
while i <= length(args)
    a = args[i]
    if a == "--seed"; global seed = parse(Int, args[i+1]); global i += 2
    elseif a == "--reps"; global reps = parse(Int, args[i+1]); global i += 2
    else; global dir = a; global i += 1
    end
end

# Pool every row (per testName) from all files of a family, keyed by nv.
# Returns Dict{nv => Vector of NamedTuple(name, tot::Dict{solver=>Float64}, its::Dict)}
function pool(fam)
    pat = Regex("^" * fam * "\\..*seed$(seed)\\.reps$(reps)\\.jld2\$")
    bysize = Dict{Int,Vector{Any}}()
    for f in sort(readdir(dir; join = true))
        occursin(pat, basename(f)) || continue
        dic = try JLD2.load(f)["dic"] catch; continue end
        haskey(dic, "testName") || continue
        names = get(dic, "names", SOLVERS)
        tns = dic["testName"]
        for (r, tn) in enumerate(tns)
            nv = haskey(dic, "nv") ? Int(dic["nv"][r]) : -1
            tot = Dict{String,Float64}(); its = Dict{String,Float64}()
            for s in SOLVERS
                s in names || continue
                tot[s] = Float64(dic["$(s)_tot"][r])
                its[s] = Float64(dic["$(s)_its"][r])
            end
            push!(get!(bysize, nv, Any[]), (name = String(tn), tot = tot, its = its))
        end
    end
    return bysize
end

ratios(rows, base, s) = [rows[i].tot[base] / rows[i].tot[s]
                         for i in eachindex(rows)
                         if haskey(rows[i].tot, base) && haskey(rows[i].tot, s) &&
                            isfinite(rows[i].tot[base]) && isfinite(rows[i].tot[s]) &&
                            rows[i].tot[s] > 0]

for fam in CHIMERA
    bysize = pool(fam)
    isempty(bysize) && (println("\n$(fam): no seed$(seed).reps$(reps) files found in $(dir)"); continue)
    println("\n=== $(fam) ===")
    for nv in sort(collect(keys(bysize)))
        rows = bysize[nv]
        n = length(rows)
        @printf("  n=%-9d  %3d samples\n", nv, n)
        for s in ["cmg-v", "cmg-k", "cmg-k-elim"]
            rs = ratios(rows, "ac", s)
            fails = count(r -> haskey(r.tot, s) && !isfinite(r.tot[s]), rows)
            if isempty(rs)
                @printf("    %-11s  fails=%d  (no converged samples)\n", s, fails)
            else
                @printf("    %-11s  fails=%-2d  ac/%s tot: min=%.3fx  med=%.3fx  max=%.3fx\n",
                        s, fails, s, minimum(rs), median(rs), maximum(rs))
            end
            # For each CMG non-convergence, did ac converge on the same graph?
            for r in rows
                if haskey(r.tot, s) && !isfinite(r.tot[s])
                    acok = haskey(r.tot, "ac") && isfinite(r.tot["ac"])
                    @printf("        NONCONV %-28s  its=%s  ac_converged=%s\n",
                            r.name, get(r.its, s, NaN), acok)
                end
            end
        end
    end
end
