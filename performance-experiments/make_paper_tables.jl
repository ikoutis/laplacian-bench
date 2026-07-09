#==========================================================
Public paper-comparison tables for the chol-vs-kcycle benchmark.

    julia --project=.. make_paper_tables.jl [DIR_OR_FILES...] \
        [--solvers ac,ac-s2m2,cmg-v,cmg-k-elim] \
        [--csv PATH] [--md PATH] [--coverage PATH]
    julia --project=.. make_paper_tables.jl --selftest

Reads every .jld2 produced by chol_vs_kcycle.jl (default directory
../performance-analyses/chol-vs-kcycle), aggregates per (file, testName) across
repetitions, and emits, for the chosen solver columns (default the paper's two
ApproxChol solvers plus CMG-legacy and CMG-K-elim):

  * paper_comparison.csv  — machine-readable, RFC-4180 quoted (testNames and
    some solver names contain commas): one row per instance, per solver the
    build/solve/tot/its/err medians and the failure count.
  * paper_comparison.md   — one table per family (median total seconds with
    median iterations in parentheses), plus a per-family accuracy line and a
    median cmg-k-elim-vs-ac total-time speedup footer. Analogous to the paper's
    per-family tables.
  * coverage.txt          — expected vs produced instance counts per family and
    per-solver failures, so "nothing was missed" is auditable.

Aggregation: median across reps of finite build/solve/tot/its; max of finite
err. Failed runs (Inf from the harness) are excluded from medians and counted.
===========================================================#

using JLD2
using Statistics
using Printf

const DEFAULT_DIR = normpath(joinpath(@__DIR__, "..", "performance-analyses", "chol-vs-kcycle"))
const DEFAULT_SOLVERS = ["ac", "ac-s2m2", "cmg-v", "cmg-k-elim"]

# Stable family order (matches benchFamilies.jl FAMILY_ORDER) for output.
const FAMILY_ORDER = ["uniform_grid", "aniso", "wgrid", "checkered",
    "sachdeva_star", "suitesparse", "spe", "chimeraIPM", "spielmanIPM",
    "uni_chimera", "uni_bndry_chimera", "wted_chimera", "wted_bndry_chimera"]

# Paper-scale instance counts for the families with a fixed sweep. IPM and SPE
# are data-dependent (full ipmMat.zip sweep / manual SPE), reported as produced.
const EXPECTED = Dict("uniform_grid" => 3, "aniso" => 7, "wgrid" => 7,
    "checkered" => 7, "sachdeva_star" => 15, "suitesparse" => 28,
    "uni_chimera" => 4, "uni_bndry_chimera" => 4, "wted_chimera" => 4,
    "wted_bndry_chimera" => 4)

fin(v) = [x for x in v if isfinite(x)]
medf(v) = isempty(v) ? Inf : median(v)
maxf(v) = isempty(v) ? Inf : maximum(v)
famkey(s) = (i = findfirst(==(s), FAMILY_ORDER); i === nothing ? length(FAMILY_ORDER) + 1 : i)
natkey(s::AbstractString) = ([parse(Int, m.match) for m in eachmatch(r"\d+", s)], String(s))

fmt(x) = !isfinite(x) ? "--" : x >= 100 ? @sprintf("%.0f", x) :
         x >= 1 ? @sprintf("%.2f", x) : @sprintf("%.3g", x)
csvnum(x) = isfinite(x) ? @sprintf("%.6g", x) : ""
function csvfield(x)
    s = x isa AbstractString ? String(x) : string(x)
    (occursin(',', s) || occursin('"', s) || occursin('\n', s)) ?
        "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end
csvrow(io, xs) = println(io, join(csvfield.(xs), ","))

# One aggregated instance: family, matrix, nv, ne, and per-solver stats.
function aggregate(files, solvers)
    recs = Any[]
    for f in files
        stored = JLD2.load(f)
        haskey(stored, "dic") || continue
        dic = stored["dic"]
        (haskey(dic, "names") && haskey(dic, "testName")) || continue
        fam = String(split(basename(f), ".")[1])
        tns = dic["testName"]
        for tn in unique(tns)
            idx = findall(==(tn), tns)
            nv = haskey(dic, "nv") ? dic["nv"][idx[1]] : -1
            ne = haskey(dic, "ne") ? dic["ne"][idx[1]] : -1
            slv = Dict{String,Any}()
            for s in solvers
                if s in dic["names"]
                    tot = dic["$(s)_tot"][idx]
                    good = fin(tot)
                    slv[s] = (build = medf(fin(dic["$(s)_build"][idx])),
                              solve = medf(fin(dic["$(s)_solve"][idx])),
                              tot = medf(good),
                              its = medf(fin(dic["$(s)_its"][idx])),
                              err = maxf(fin(dic["$(s)_err"][idx])),
                              fail = length(tot) - length(good),
                              runs = length(idx),
                              ran = true)
                else
                    slv[s] = (build = NaN, solve = NaN, tot = NaN, its = NaN,
                              err = NaN, fail = 0, runs = 0, ran = false)
                end
            end
            push!(recs, (family = fam, matrix = String(tn), nv = nv, ne = ne, slv = slv))
        end
    end
    sort!(recs; by = r -> (famkey(r.family), natkey(r.matrix)))
    return recs
end

function write_csv(path, recs, solvers)
    open(path, "w") do io
        header = ["family", "matrix", "nv", "ne"]
        for s in solvers
            append!(header, ["$(s)_build_s", "$(s)_solve_s", "$(s)_tot_s",
                             "$(s)_its", "$(s)_err", "$(s)_fail"])
        end
        csvrow(io, header)
        for r in recs
            row = Any[r.family, r.matrix, r.nv, r.ne]
            for s in solvers
                t = r.slv[s]
                append!(row, [csvnum(t.build), csvnum(t.solve), csvnum(t.tot),
                              csvnum(t.its), csvnum(t.err), t.fail])
            end
            csvrow(io, row)
        end
    end
end

function write_md(path, recs, solvers)
    open(path, "w") do io
        println(io, "# CMG vs ApproxChol — benchmark comparison\n")
        println(io, "Median total time in seconds with median iteration count in ",
                "parentheses; `--` = solver not run, `FAIL` = all reps failed, ",
                "`[Nf]` = N of the reps failed. One table per benchmark family.\n")
        fams = unique([r.family for r in recs])
        sort!(fams; by = famkey)
        for fam in fams
            frecs = [r for r in recs if r.family == fam]
            println(io, "## $(fam)  ($(length(frecs)) instances)\n")
            header = vcat(["matrix", "nv", "ne"], ["$(s)" for s in solvers])
            println(io, "| ", join(header, " | "), " |")
            println(io, "|", repeat("---|", length(header)))
            for r in frecs
                cells = String[]
                for s in solvers
                    t = r.slv[s]
                    if !t.ran
                        push!(cells, "--")
                    elseif !isfinite(t.tot)
                        push!(cells, "FAIL")
                    else
                        push!(cells, string(fmt(t.tot), " (", fmt(t.its), ")",
                            t.fail > 0 ? " [$(t.fail)f]" : ""))
                    end
                end
                println(io, "| ", join(vcat([r.matrix, string(r.nv), string(r.ne)], cells), " | "), " |")
            end
            # per-family accuracy (max relres) and cmg-k-elim vs ac speedup
            accs = [string(s, " ", fmt(maxf([r.slv[s].err for r in frecs if isfinite(r.slv[s].err)]))) for s in solvers]
            println(io, "\nMax relres: ", join(accs, ", "), ".")
            if "ac" in solvers && "cmg-k-elim" in solvers
                sp = fin([r.slv["ac"].tot / r.slv["cmg-k-elim"].tot for r in frecs
                          if isfinite(r.slv["ac"].tot) && isfinite(r.slv["cmg-k-elim"].tot) && r.slv["cmg-k-elim"].tot > 0])
                if !isempty(sp)
                    println(io, "Median cmg-k-elim vs ac total-time speedup: ",
                            @sprintf("%.2fx", median(sp)), " over $(length(sp)) instances.")
                end
            end
            println(io)
        end
    end
end

function write_coverage(path, recs, solvers)
    open(path, "w") do io
        produced = Dict{String,Int}()
        for r in recs
            produced[r.family] = get(produced, r.family, 0) + 1
        end
        println(io, "Coverage report — expected vs produced instances per family\n")
        for fam in FAMILY_ORDER
            got = get(produced, fam, 0)
            if haskey(EXPECTED, fam)
                exp = EXPECTED[fam]
                mark = got == exp ? "OK" : got == 0 ? "MISSING" : "PARTIAL"
                println(io, @sprintf("  %-20s %3d / %-3d  %s", fam, got, exp, mark))
            elseif fam == "spe"
                println(io, @sprintf("  %-20s %3d / manual  %s", fam, got,
                    got == 0 ? "skipped (SPE needs manual spe.zip)" : "present"))
            else  # chimeraIPM / spielmanIPM — full sweep is data-dependent
                println(io, @sprintf("  %-20s %3d / (sweep)  %s", fam, got,
                    got == 0 ? "MISSING (stage ipmMat.zip or run download_data.jl)" : "present"))
            end
        end
        println(io, "\nPer-solver failed instances (all reps Inf):")
        for s in solvers
            fams_failed = String[]
            for fam in FAMILY_ORDER
                nf = count(r -> r.family == fam && r.slv[s].ran && !isfinite(r.slv[s].tot), recs)
                nf > 0 && push!(fams_failed, "$(fam):$(nf)")
            end
            println(io, @sprintf("  %-12s %s", s, isempty(fams_failed) ? "none" : join(fams_failed, " ")))
        end
    end
end

# ---------------------------------------------------------------- self-test
function selftest()
    mktempdir() do dir
        dic = Dict{String,Any}()
        dic["names"] = ["ac", "cmg-k-elim"]
        # two instances (one testName contains a comma), 3 reps each
        dic["testName"] = ["uc.i1,eps0.1", "uc.i1,eps0.1", "uc.i1,eps0.1", "gridB", "gridB", "gridB"]
        dic["nv"] = [10, 10, 10, 20, 20, 20]
        dic["ne"] = [30, 30, 30, 60, 60, 60]
        dic["ac_build"] = [1.0, 1.0, 1.0, 2.0, 2.0, 2.0]
        dic["ac_solve"] = [0.5, 0.5, 0.5, 1.0, 1.0, 1.0]
        dic["ac_tot"] = [1.5, 1.5, 1.5, 3.0, 3.0, 3.0]
        dic["ac_its"] = [7.0, 7.0, 7.0, 8.0, 8.0, 8.0]
        dic["ac_err"] = [1e-9, 1e-9, 1e-9, 2e-9, 2e-9, 2e-9]
        # cmg-k-elim: 3x faster on inst 1; all-fail on inst 2 (Inf)
        dic["cmg-k-elim_build"] = [0.3, 0.3, 0.3, Inf, Inf, Inf]
        dic["cmg-k-elim_solve"] = [0.2, 0.2, 0.2, Inf, Inf, Inf]
        dic["cmg-k-elim_tot"] = [0.5, 0.5, 0.5, Inf, Inf, Inf]
        dic["cmg-k-elim_its"] = [1.0, 1.0, 1.0, Inf, Inf, Inf]
        dic["cmg-k-elim_err"] = [1e-11, 1e-11, 1e-11, Inf, Inf, Inf]
        f = joinpath(dir, "chimeraIPM.paper.seed1.reps3.jld2")
        JLD2.save(f, "dic", dic)

        solvers = ["ac", "cmg-k-elim"]
        recs = aggregate([f], solvers)
        @assert length(recs) == 2 "expected 2 instances, got $(length(recs))"
        r1 = recs[findfirst(r -> r.matrix == "uc.i1,eps0.1", recs)]
        @assert r1.slv["ac"].tot == 1.5
        @assert r1.slv["cmg-k-elim"].tot == 0.5
        @assert r1.slv["cmg-k-elim"].fail == 0
        r2 = recs[findfirst(r -> r.matrix == "gridB", recs)]
        @assert !isfinite(r2.slv["cmg-k-elim"].tot) "all-fail must aggregate to Inf"
        @assert r2.slv["cmg-k-elim"].fail == 3

        csv = joinpath(dir, "t.csv"); write_csv(csv, recs, solvers)
        md = joinpath(dir, "t.md"); write_md(md, recs, solvers)
        cov = joinpath(dir, "cov.txt"); write_coverage(cov, recs, solvers)
        csvtxt = read(csv, String)
        @assert occursin("\"uc.i1,eps0.1\"", csvtxt) "comma testName must be quoted in CSV"
        @assert occursin("chimeraIPM,\"uc.i1,eps0.1\",10,30,1,0.5,1.5,7,1e-09,0", csvtxt) "ac row values"
        mdtxt = read(md, String)
        @assert occursin("FAIL", mdtxt) "all-fail instance must show FAIL in MD"
        @assert occursin("speedup: 3.00x", mdtxt) "speedup footer (1.5/0.5 = 3x)"
        println("make_paper_tables selftest: PASSED")
    end
end

# ------------------------------------------------------------------- main
function main()
    if "--selftest" in ARGS
        selftest()
        return
    end
    inputs = String[]
    solvers = copy(DEFAULT_SOLVERS)
    csvpath = nothing; mdpath = nothing; covpath = nothing
    let i = 1
        while i <= length(ARGS)
            a = ARGS[i]
            if a == "--solvers" && i < length(ARGS)
                solvers = String.(split(ARGS[i+1], ',')); i += 2
            elseif a == "--csv" && i < length(ARGS)
                csvpath = ARGS[i+1]; i += 2
            elseif a == "--md" && i < length(ARGS)
                mdpath = ARGS[i+1]; i += 2
            elseif a == "--coverage" && i < length(ARGS)
                covpath = ARGS[i+1]; i += 2
            else
                push!(inputs, a); i += 1
            end
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
        return
    end
    outdir = dirname(files[1])
    csvpath = csvpath === nothing ? joinpath(outdir, "paper_comparison.csv") : csvpath
    mdpath = mdpath === nothing ? joinpath(outdir, "paper_comparison.md") : mdpath
    covpath = covpath === nothing ? joinpath(outdir, "coverage.txt") : covpath

    recs = aggregate(files, solvers)
    write_csv(csvpath, recs, solvers)
    write_md(mdpath, recs, solvers)
    write_coverage(covpath, recs, solvers)
    println("solvers: ", join(solvers, ", "))
    println("instances: ", length(recs), " across ", length(files), " result files")
    println("wrote ", csvpath)
    println("wrote ", mdpath)
    println("wrote ", covpath)
    println("\n--- coverage ---")
    print(read(covpath, String))
end

main()
