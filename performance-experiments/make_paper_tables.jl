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
    some solver names contain commas): one row per instance (all instances,
    unfiltered), with the system kind (lap/sddm), per solver the
    build/solve/tot/its/err medians and the failure count.
  * paper_comparison.md   — a leading speedup-summary section (per-family
    median total AND solve-only speedups of cmg-k-elim vs ac and cmg-v, plus a
    worst-case-per-family table naming the worst instance), then one table per
    family (median total seconds with median iterations in parentheses), a
    per-family accuracy line, and a median cmg-k-elim-vs-ac total-time speedup
    footer. SuiteSparse is split into a Laplacian and an SDDM table with
    ne > 1000, mirroring the paper.
  * coverage.txt          — expected vs produced instance counts per family
    (SuiteSparse also shows the Laplacian/SDDM split) and per-solver failures,
    so "nothing was missed" is auditable.

Aggregation: median across reps of finite build/solve/tot/its and err (matching
the paper's printMedian). Failed runs (Inf from the harness) are excluded from
the medians and counted separately. The chimera families draw a fresh random
graph each rep, so their per-size samples are collapsed into one median row per
size (matching the paper's per-size chimera tables).
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

# The paper's suitesparse.ipynb keeps only matrices with ne > 1000 in its tables
# and splits them into a Laplacian and an SDDM table (partitioning on whether
# the diagonal is zero). We mirror both for the `suitesparse` family.
const NE_TABLE_MIN = 1000

# Chimera families draw a NEW random graph each rep (the instance index derives
# from baseseed+rep), so a size run yields several distinct samples that all
# share nv (= the size). We collapse those samples per size into one median row,
# matching the paper's per-size chimera tables.
const CHIMERA_FAMILIES = ["uni_chimera", "uni_bndry_chimera",
    "wted_chimera", "wted_bndry_chimera"]

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

# Collapse a group of chimera samples (same family + nv) into one rec: median
# over the samples of each per-solver stat, median ne, summed failures.
function merge_samples(group, solvers)
    fam = group[1].family
    nv = group[1].nv
    ne = round(Int, medf([Float64(r.ne) for r in group]))
    n = length(group)
    slv = Dict{String,Any}()
    for s in solvers
        ts = [r.slv[s] for r in group if r.slv[s].ran]
        if isempty(ts)
            slv[s] = (build = NaN, solve = NaN, tot = NaN, its = NaN,
                      err = NaN, fail = 0, runs = 0, ran = false)
        else
            slv[s] = (build = medf(fin([t.build for t in ts])),
                      solve = medf(fin([t.solve for t in ts])),
                      tot   = medf(fin([t.tot for t in ts])),
                      its   = medf(fin([t.its for t in ts])),
                      err   = medf(fin([t.err for t in ts])),
                      fail  = sum(t.fail for t in ts),
                      runs  = sum(t.runs for t in ts),
                      ran   = true)
        end
    end
    return (family = fam, matrix = "$(fam)(n=$(nv), $(n) samples)",
            nv = nv, ne = ne, kind = "lap", slv = slv)
end

# Collapse the chimera families' per-sample recs into one row per (family, nv).
function collapse_chimera(recs, solvers)
    out = Any[]
    groups = Dict{Tuple{String,Int},Vector{Any}}()
    order = Tuple{String,Int}[]
    for r in recs
        if r.family in CHIMERA_FAMILIES
            key = (r.family, r.nv)
            haskey(groups, key) || (groups[key] = Any[]; push!(order, key))
            push!(groups[key], r)
        else
            push!(out, r)
        end
    end
    for key in order
        push!(out, merge_samples(groups[key], solvers))
    end
    return out
end

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
            # System kind (:lap/:sddm), recorded per row by chol_vs_kcycle.jl.
            # Absent in older result files -> "" (no SuiteSparse Lap/SDDM split).
            kind = haskey(dic, "kind") ? String(dic["kind"][idx[1]]) : ""
            slv = Dict{String,Any}()
            for s in solvers
                if s in dic["names"]
                    tot = dic["$(s)_tot"][idx]
                    good = fin(tot)
                    slv[s] = (build = medf(fin(dic["$(s)_build"][idx])),
                              solve = medf(fin(dic["$(s)_solve"][idx])),
                              tot = medf(good),
                              its = medf(fin(dic["$(s)_its"][idx])),
                              err = medf(fin(dic["$(s)_err"][idx])),
                              fail = length(tot) - length(good),
                              runs = length(idx),
                              ran = true)
                else
                    slv[s] = (build = NaN, solve = NaN, tot = NaN, its = NaN,
                              err = NaN, fail = 0, runs = 0, ran = false)
                end
            end
            push!(recs, (family = fam, matrix = String(tn), nv = nv, ne = ne, kind = kind, slv = slv))
        end
    end
    recs = collapse_chimera(recs, solvers)
    sort!(recs; by = r -> (famkey(r.family), natkey(r.matrix)))
    return recs
end

function write_csv(path, recs, solvers)
    open(path, "w") do io
        header = ["family", "matrix", "nv", "ne", "kind"]
        for s in solvers
            append!(header, ["$(s)_build_s", "$(s)_solve_s", "$(s)_tot_s",
                             "$(s)_its", "$(s)_err", "$(s)_fail"])
        end
        csvrow(io, header)
        for r in recs
            row = Any[r.family, r.matrix, r.nv, r.ne, r.kind]
            for s in solvers
                t = r.slv[s]
                append!(row, [csvnum(t.build), csvnum(t.solve), csvnum(t.tot),
                              csvnum(t.its), csvnum(t.err), t.fail])
            end
            csvrow(io, row)
        end
    end
end

# ------------------------------------------------ speedup summary section

const SPEEDUP_TARGET = "cmg-k-elim"
# cmg-v isolates the cycle change; cmg-k isolates the elimination (k vs k-elim).
# Bases absent from a run are simply omitted from the summary.
const SPEEDUP_BASES = ["ac", "cmg-v", "cmg-k"]

# Group recs for the summary: FAMILY_ORDER order, with suitesparse split into
# its Laplacian and SDDM halves when the kind column is recorded.
function summary_groups(recs)
    groups = Tuple{String,Vector{Any}}[]
    fams = unique([r.family for r in recs])
    sort!(fams; by = famkey)
    for fam in fams
        frecs = [r for r in recs if r.family == fam]
        if fam == "suitesparse" && any(!isempty(r.kind) for r in frecs)
            lap = [r for r in frecs if r.kind == "lap"]
            sddm = [r for r in frecs if r.kind == "sddm"]
            isempty(lap) || push!(groups, ("suitesparse (Laplacian)", lap))
            isempty(sddm) || push!(groups, ("suitesparse (SDDM)", sddm))
        else
            push!(groups, (fam, frecs))
        end
    end
    return groups
end

# Per-instance speedup base/target on `field` (:tot or :solve), with labels.
function speedup_pairs(rows, base, field)
    out = Tuple{Float64,String}[]
    for r in rows
        (r.slv[base].ran && r.slv[SPEEDUP_TARGET].ran) || continue
        x = getfield(r.slv[base], field)
        y = getfield(r.slv[SPEEDUP_TARGET], field)
        (isfinite(x) && isfinite(y) && y > 0) || continue
        push!(out, (x / y, r.matrix))
    end
    return out
end

spx(x) = string(fmt(x), "x")

# Median (total + solve) and worst-case per-family speedup tables for
# SPEEDUP_TARGET against whichever of SPEEDUP_BASES are in the run. No-op when
# the target or every baseline is absent.
function write_speedup_summary(io, recs, solvers)
    SPEEDUP_TARGET in solvers || return
    bases = [b for b in SPEEDUP_BASES if b in solvers]
    isempty(bases) && return
    groups = summary_groups(recs)

    println(io, "## Speedup summary — `$(SPEEDUP_TARGET)` vs ",
            join(("`$(b)`" for b in bases), ", "), "\n")
    println(io, "Median per-family speedup of `$(SPEEDUP_TARGET)` (>1 = ",
            "`$(SPEEDUP_TARGET)` faster). `total` = build + solve; `solve` ",
            "excludes the one-time build (the relevant number when one ",
            "factorization serves many right-hand sides).\n")
    header = vcat(["family", "instances"],
                  reduce(vcat, [["vs $(b) total", "vs $(b) solve"] for b in bases]))
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", repeat("---|", length(header)))
    for (label, rows) in groups
        cells = String[label, string(length(rows))]
        for b in bases
            t = speedup_pairs(rows, b, :tot)
            s = speedup_pairs(rows, b, :solve)
            push!(cells, isempty(t) ? "--" : spx(median(first.(t))))
            push!(cells, isempty(s) ? "--" : spx(median(first.(s))))
        end
        println(io, "| ", join(cells, " | "), " |")
    end

    println(io, "\nWorst case per family — minimum per-instance total-time ",
            "speedup, and the instance it occurs on:\n")
    header = vcat(["family"], reduce(vcat, [["vs $(b) worst", "instance"] for b in bases]))
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", repeat("---|", length(header)))
    for (label, rows) in groups
        cells = String[label]
        for b in bases
            t = speedup_pairs(rows, b, :tot)
            if isempty(t)
                append!(cells, ["--", "--"])
            else
                v, m = minimum(t)
                append!(cells, [spx(v), m])
            end
        end
        println(io, "| ", join(cells, " | "), " |")
    end
    println(io)
end

# Emit one Markdown table (header + rows) plus a per-table accuracy line and a
# median cmg-k-elim-vs-ac total-time speedup footer, for the given rows.
function emit_md_table(io, title, rows, solvers)
    println(io, "## $(title)  ($(length(rows)) instances)\n")
    header = vcat(["matrix", "nv", "ne"], ["$(s)" for s in solvers])
    println(io, "| ", join(header, " | "), " |")
    println(io, "|", repeat("---|", length(header)))
    for r in rows
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
    isempty(rows) && (println(io); return)
    # per-table accuracy (max relres) and cmg-k-elim vs ac speedup
    accs = [string(s, " ", fmt(maxf([r.slv[s].err for r in rows if isfinite(r.slv[s].err)]))) for s in solvers]
    println(io, "\nMax relres: ", join(accs, ", "), ".")
    if "ac" in solvers && "cmg-k-elim" in solvers
        sp = fin([r.slv["ac"].tot / r.slv["cmg-k-elim"].tot for r in rows
                  if isfinite(r.slv["ac"].tot) && isfinite(r.slv["cmg-k-elim"].tot) && r.slv["cmg-k-elim"].tot > 0])
        if !isempty(sp)
            println(io, "Median cmg-k-elim vs ac total-time speedup: ",
                    @sprintf("%.2fx", median(sp)), " over $(length(sp)) instances.")
        end
    end
    println(io)
end

function write_md(path, recs, solvers)
    open(path, "w") do io
        println(io, "# CMG vs ApproxChol — benchmark comparison\n")
        println(io, "Median total time in seconds with median iteration count in ",
                "parentheses; `--` = solver not run, `FAIL` = all reps failed, ",
                "`[Nf]` = N of the reps failed. One table per benchmark family; ",
                "SuiteSparse is split into Laplacian and SDDM tables (ne > ",
                "$(NE_TABLE_MIN)) as in the paper.\n")
        write_speedup_summary(io, recs, solvers)
        fams = unique([r.family for r in recs])
        sort!(fams; by = famkey)
        for fam in fams
            frecs = [r for r in recs if r.family == fam]
            if fam == "suitesparse" && any(!isempty(r.kind) for r in frecs)
                big = [r for r in frecs if r.ne > NE_TABLE_MIN]
                small = [r for r in frecs if r.ne <= NE_TABLE_MIN]
                lap = [r for r in big if r.kind == "lap"]
                sddm = [r for r in big if r.kind == "sddm"]
                other = [r for r in big if r.kind != "lap" && r.kind != "sddm"]
                emit_md_table(io, "suitesparse — Laplacian (ne > $(NE_TABLE_MIN))", lap, solvers)
                emit_md_table(io, "suitesparse — SDDM (ne > $(NE_TABLE_MIN))", sddm, solvers)
                isempty(other) || emit_md_table(io, "suitesparse — unlabeled kind (ne > $(NE_TABLE_MIN))", other, solvers)
                if !isempty(small)
                    println(io, "_Excluded from the SuiteSparse tables (ne ≤ $(NE_TABLE_MIN), as in the paper): ",
                            join([string(r.matrix, " (ne=", r.ne, ")") for r in small], ", "), "._\n")
                end
            else
                emit_md_table(io, fam, frecs, solvers)
            end
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
                # SuiteSparse: show the Laplacian/SDDM split and the ne<=1000
                # matrices the paper drops from its tables (still counted here).
                if fam == "suitesparse"
                    frecs = [r for r in recs if r.family == "suitesparse"]
                    if any(!isempty(r.kind) for r in frecs)
                        nlap = count(r -> r.kind == "lap", frecs)
                        nsddm = count(r -> r.kind == "sddm", frecs)
                        nsmall = count(r -> r.ne <= NE_TABLE_MIN, frecs)
                        println(io, @sprintf("  %-20s   Laplacian %d, SDDM %d; %d with ne<=%d excluded from tables",
                            "", nlap, nsddm, nsmall, NE_TABLE_MIN))
                    end
                end
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
        # no "kind" recorded in this dic -> empty kind field between ne and the solver cols
        @assert occursin("chimeraIPM,\"uc.i1,eps0.1\",10,30,,1,0.5,1.5,7,1e-09,0", csvtxt) "ac row values"
        mdtxt = read(md, String)
        @assert occursin("FAIL", mdtxt) "all-fail instance must show FAIL in MD"
        @assert occursin("speedup: 3.00x", mdtxt) "speedup footer (1.5/0.5 = 3x)"
        # speedup summary: vs-ac only (no cmg-v in solvers); the all-fail
        # instance is excluded, so median == worst == the single valid instance
        @assert occursin("Speedup summary", mdtxt) "summary section present"
        @assert occursin("vs ac total", mdtxt) && !occursin("vs cmg-v", mdtxt) &&
                !occursin("vs cmg-k ", mdtxt) "only ran baselines appear"
        @assert occursin("| chimeraIPM | 2 | 3.00x | 2.50x |", mdtxt) "median total 1.5/0.5, solve 0.5/0.2"
        @assert occursin("| chimeraIPM | 3.00x | uc.i1,eps0.1 |", mdtxt) "worst-case row names the instance"

        # --- SuiteSparse Lap/SDDM split + ne>1000 filter ---
        ss = Dict{String,Any}()
        ss["names"] = ["ac", "cmg-k-elim"]
        ss["testName"] = ["McRae/ecology1", "HB/bcsstm08", "HB/nos6"]
        ss["nv"] = [1000, 900, 675]
        ss["ne"] = [2000, 5000, 500]          # nos6 below threshold -> excluded from tables
        ss["kind"] = ["lap", "sddm", "sddm"]
        for s in ("ac", "cmg-k-elim"), col in ("build", "solve", "tot", "its", "err")
            ss["$(s)_$(col)"] = [1.0, 1.0, 1.0]
        end
        fss = joinpath(dir, "suitesparse.paper.seed1.reps1.jld2")
        JLD2.save(fss, "dic", ss)
        ssrecs = aggregate([fss], solvers)
        @assert length(ssrecs) == 3 "expected 3 suitesparse instances"
        @assert ssrecs[findfirst(r -> r.matrix == "McRae/ecology1", ssrecs)].kind == "lap" "kind captured"
        sscsv = joinpath(dir, "ss.csv"); write_csv(sscsv, ssrecs, solvers)
        ssmd  = joinpath(dir, "ss.md");  write_md(ssmd, ssrecs, solvers)
        sscov = joinpath(dir, "ss.cov"); write_coverage(sscov, ssrecs, solvers)
        sscsvt = read(sscsv, String); ssmdt = read(ssmd, String); sscovt = read(sscov, String)
        @assert occursin("suitesparse,McRae/ecology1,1000,2000,lap,", sscsvt) "kind column in CSV"
        @assert occursin("suitesparse — Laplacian", ssmdt) "Laplacian subtable header"
        @assert occursin("suitesparse — SDDM", ssmdt) "SDDM subtable header"
        @assert occursin("| McRae/ecology1 |", ssmdt) "lap ne>1000 row present"
        @assert occursin("| HB/bcsstm08 |", ssmdt) "sddm ne>1000 row present"
        @assert !occursin("| HB/nos6 |", ssmdt) "ne<=1000 excluded from tables"
        @assert occursin("HB/nos6 (ne=500)", ssmdt) "excluded note lists the small matrix"
        @assert occursin("Laplacian 1, SDDM 2", sscovt) "coverage shows Lap/SDDM split"

        # --- chimera per-size collapse (random samples -> one row per size) ---
        ch = Dict{String,Any}()
        ch["names"] = ["ac", "cmg-k-elim"]
        # two sizes; size 10000 has two random samples, size 20000 has one
        ch["testName"] = ["uni_chimera(10000 1001)", "uni_chimera(10000 1002)", "uni_chimera(20000 1001)"]
        ch["nv"] = [10000, 10000, 20000]
        ch["ne"] = [100, 200, 300]
        ch["kind"] = ["lap", "lap", "lap"]
        for s in ("ac", "cmg-k-elim"), col in ("build", "solve", "its", "err")
            ch["$(s)_$(col)"] = [1.0, 1.0, 1.0]
        end
        ch["ac_tot"] = [1.0, 3.0, 9.0]         # size 10000 -> median(1,3)=2
        ch["cmg-k-elim_tot"] = [1.0, 1.0, 1.0]
        fch = joinpath(dir, "uni_chimera.n1e4.seed1.reps3.jld2")
        JLD2.save(fch, "dic", ch)
        chrecs = aggregate([fch], solvers)
        @assert length(chrecs) == 2 "chimera should collapse 3 samples -> 2 size rows, got $(length(chrecs))"
        c10 = chrecs[findfirst(r -> r.nv == 10000, chrecs)]
        @assert c10.slv["ac"].tot == 2.0 "size 10000 tot = median(1,3) = 2"
        @assert c10.ne == 150 "size 10000 ne = median(100,200) = 150"
        @assert occursin("2 samples", c10.matrix) "collapsed label notes sample count"
        @assert count(r -> r.family == "uni_chimera", chrecs) == 2 "2 collapsed uni_chimera size rows"
        chcov = joinpath(dir, "ch.cov"); write_coverage(chcov, chrecs, solvers)
        @assert occursin(r"uni_chimera\s+2 / 4", read(chcov, String)) "coverage counts collapsed size rows (2/4)"

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
