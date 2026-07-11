#!/usr/bin/env julia
# Why is CMG sometimes much slower than ac on a converged chimera draw?
#
# The tail analysis showed converged draws where ac_tot/cmg-k-elim_tot is very
# small (ac 15-50x faster). This finds, per (family, size), the worst such
# *converged* instances and prints the build/solve/iteration breakdown for ac,
# cmg-k, and cmg-k-elim — so we can see whether CMG's slowness is build-time
# (hierarchy setup) or solve-time (many iterations / high condition number).
#
# Usage (compute node, after `source wulver/env_wulver.sh`):
#   julia --project=.. analyze_chimera_slow.jl [RESULTS_DIR] [--seed S] [--reps R] [--top K]
# Defaults: RESULTS_DIR=../performance-analyses/chol-vs-kcycle, seed 2, reps 1, top 3.

using JLD2, Printf

const CHIMERA = ["uni_chimera", "uni_bndry_chimera", "wted_chimera", "wted_bndry_chimera"]
const SOLVERS = ["ac", "cmg-k", "cmg-k-elim"]
const FIELDS  = ["build", "solve", "tot", "its"]

dir = "../performance-analyses/chol-vs-kcycle"
seed, reps, topk = 2, 1, 3
i = 1
while i <= length(ARGS)
    a = ARGS[i]
    if a == "--seed"; global seed = parse(Int, ARGS[i+1]); global i += 2
    elseif a == "--reps"; global reps = parse(Int, ARGS[i+1]); global i += 2
    elseif a == "--top"; global topk = parse(Int, ARGS[i+1]); global i += 2
    else; global dir = a; global i += 1
    end
end

# Pool rows (per testName) across a family's files, keyed by nv. Each row carries
# per-solver (build, solve, tot, its).
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
ratio(row, s) = conv(row, "ac") && conv(row, s) ? row.s["ac"]["tot"] / row.s[s]["tot"] : Inf

function line(row, s)
    haskey(row.s, s) || return @sprintf("    %-11s (not run)", s)
    d = row.s[s]
    bshare = isfinite(d["tot"]) && d["tot"] > 0 ? 100 * d["build"] / d["tot"] : NaN
    conv_mark = isfinite(d["tot"]) ? "" : "  [DID NOT CONVERGE]"
    @printf("    %-11s build=%.3g  solve=%.3g  tot=%.3g  its=%.0f  (build %.0f%%)%s",
            s, d["build"], d["solve"], d["tot"], d["its"], bshare, conv_mark)
    return ""
end

for fam in CHIMERA
    bysize = pool(fam)
    isempty(bysize) && (println("\n$(fam): no seed$(seed).reps$(reps) files in $(dir)"); continue)
    println("\n=== $(fam) ===")
    for nv in sort(collect(keys(bysize)))
        rows = bysize[nv]
        # converged-in-cmg-k-elim rows, worst (smallest ac/elim ratio) first
        cr = filter(r -> conv(r, "cmg-k-elim"), rows)
        sort!(cr; by = r -> ratio(r, "cmg-k-elim"))
        @printf("  n=%-9d  %d samples (%d converged in cmg-k-elim)\n", nv, length(rows), length(cr))
        for r in cr[1:min(topk, end)]
            @printf("   %-26s  ac/cmg-k-elim tot = %.3fx\n", r.name, ratio(r, "cmg-k-elim"))
            for s in SOLVERS
                msg = line(r, s); isempty(msg) ? println() : println(msg)
            end
        end
    end
end
