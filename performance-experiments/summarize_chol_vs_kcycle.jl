#==========================================================
Summarize chol-vs-kcycle results.

    julia --project=.. summarize_chol_vs_kcycle.jl [DIR_OR_FILES...] [--tsv PATH]

Reads every .jld2 produced by chol_vs_kcycle.jl (default directory:
../performance-analyses/chol-vs-kcycle), aggregates per (file, testName,
solver) across repetitions, prints per-file Markdown tables of median total
time (median iterations in parentheses), and writes one machine-readable
TSV (tab-separated, since testNames contain commas).

TSV columns: file, testName, solver, runs, failures,
             med_build_s, med_solve_s, med_tot_s, med_its, max_err
Failed runs (recorded as Inf by the harness) are excluded from the medians
and counted in `failures`.
===========================================================#

using JLD2
using Statistics
using Printf
using DelimitedFiles

const DEFAULT_DIR = normpath(joinpath(@__DIR__, "..", "performance-analyses", "chol-vs-kcycle"))

inputs = String[]
tsvpath = nothing
let i = 1
    while i <= length(ARGS)
        if ARGS[i] == "--tsv" && i < length(ARGS)
            global tsvpath = ARGS[i+1]
            i += 1
        else
            push!(inputs, ARGS[i])
        end
        i += 1
    end
end
isempty(inputs) && push!(inputs, DEFAULT_DIR)

files = String[]
for inp in inputs
    if isdir(inp)
        append!(files, sort!([joinpath(inp, f) for f in readdir(inp) if endswith(f, ".jld2")]))
    elseif isfile(inp)
        push!(files, inp)
    else
        @warn "no such file or directory: $(inp)"
    end
end
if isempty(files)
    println("no .jld2 result files found (looked at: $(join(inputs, ", ")))")
    exit(0)
end
tsvpath = tsvpath === nothing ? joinpath(dirname(files[1]), "summary.tsv") : tsvpath

fin(v) = [x for x in v if isfinite(x)]
fmt(x) = !isfinite(x) ? "--" : x >= 100 ? @sprintf("%.0f", x) :
         x >= 1 ? @sprintf("%.2f", x) : @sprintf("%.3g", x)

rows = Vector{Any}[]

for f in files
    stored = JLD2.load(f)
    haskey(stored, "dic") || (println("skipping $(basename(f)) (no `dic` variable)"); continue)
    dic = stored["dic"]
    (haskey(dic, "names") && haskey(dic, "testName")) ||
        (println("skipping $(basename(f)) (not a results dict)"); continue)

    solvers = dic["names"]
    tns = dic["testName"]
    utns = unique(tns)

    println()
    println("## $(basename(f))")
    println()
    header = vcat(["testName", "nv"], ["$(s) tot(its)" for s in solvers])
    println("| ", join(header, " | "), " |")
    println("|", repeat("---|", length(header)))

    for tn in utns
        idx = findall(==(tn), tns)
        nv = haskey(dic, "nv") ? dic["nv"][idx[1]] : -1
        cells = String[]
        for s in solvers
            tot = dic["$(s)_tot"][idx]
            its = dic["$(s)_its"][idx]
            good = fin(tot)
            nfail = length(tot) - length(good)
            cell = isempty(good) ? "FAIL" :
                string(fmt(median(good)), " (", fmt(median(fin(its))), ")",
                       nfail > 0 ? " [$(nfail)f]" : "")
            push!(cells, cell)

            build = fin(dic["$(s)_build"][idx])
            solve = fin(dic["$(s)_solve"][idx])
            err = fin(dic["$(s)_err"][idx])
            push!(rows, Any[basename(f), tn, s, length(idx), nfail,
                isempty(build) ? Inf : median(build),
                isempty(solve) ? Inf : median(solve),
                isempty(good) ? Inf : median(good),
                isempty(fin(its)) ? Inf : median(fin(its)),
                isempty(err) ? Inf : maximum(err)])
        end
        println("| ", join(vcat([tn, string(nv)], cells), " | "), " |")
    end
end

open(tsvpath, "w") do io
    writedlm(io, [["file" "testName" "solver" "runs" "failures" "med_build_s" "med_solve_s" "med_tot_s" "med_its" "max_err"]], '\t')
    for r in rows
        writedlm(io, [permutedims(r)], '\t')
    end
end
println()
println("wrote $(tsvpath) ($(length(rows)) rows)")
