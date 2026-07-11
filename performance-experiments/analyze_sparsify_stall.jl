#==========================================================
Per-sample stalled-draw analysis for the sparsify-on-stall comparison.

The paper tables collapse each chimera family's samples into ONE median row
per size, which mixes the (majority) non-stalling draws with the (minority)
stalling ones — so the median hides exactly the draws sparsify is meant to
rescue. This script instead works per SAMPLE: it splits every draw into
"stalled" (some sparsify column injected a level, `_inj > 0`) vs "clean", and
reports the head-to-head on the stalled draws only.

    julia --project=.. analyze_sparsify_stall.jl [DIR_OR_FILES...]
        [--solvers ac,cmg-k-elim,cmg-sparsify-l,cmg-sparsify-ks]

Default dir: ../performance-analyses/chol-vs-kcycle. For each family it prints,
over the stalled draws: how many stalled, median total-time and iterations per
solver, and the per-draw median speedup of each sparsify column vs cmg-k-elim
and vs ac (>1 = sparsify faster). The clean draws are shown for contrast.
===========================================================#

using JLD2
using Statistics
using Printf

const DEFAULT_DIR = normpath(joinpath(@__DIR__, "..", "performance-analyses", "chol-vs-kcycle"))
const DEFAULT_SOLVERS = ["ac", "cmg-k-elim", "cmg-sparsify-l", "cmg-sparsify-ks"]
const SPARSIFY = ["cmg-sparsify-l", "cmg-sparsify-ks"]
const REF = "cmg-k-elim"

fin(v) = [x for x in v if isfinite(x)]
medf(v) = isempty(v) ? NaN : median(v)
fmt(x) = !isfinite(x) ? "--" : x >= 100 ? @sprintf("%.0f", x) :
         x >= 1 ? @sprintf("%.2f", x) : @sprintf("%.3g", x)

# Per-sample records: one dict per (file row), carrying family + per-solver
# tot/its/inj (Inf/absent -> skipped downstream via fin()).
function load_rows(files, solvers)
    rows = Vector{Dict{String,Any}}()
    for f in files
        stored = try
            JLD2.load(f)
        catch e
            @warn "skip unreadable file" file = f err = e
            continue
        end
        haskey(stored, "dic") || continue
        dic = stored["dic"]
        (haskey(dic, "names") && haskey(dic, "testName")) || continue
        fam = String(split(basename(f), ".")[1])
        names = dic["names"]
        n = length(dic["testName"])
        for i in 1:n
            rec = Dict{String,Any}("family" => fam,
                                   "matrix" => String(dic["testName"][i]),
                                   "nv" => haskey(dic, "nv") ? dic["nv"][i] : -1)
            for s in solvers
                if s in names
                    rec["$(s)_tot"] = Float64(dic["$(s)_tot"][i])
                    rec["$(s)_its"] = Float64(dic["$(s)_its"][i])
                    rec["$(s)_inj"] = haskey(dic, "$(s)_inj") ? Float64(dic["$(s)_inj"][i]) : 0.0
                else
                    rec["$(s)_tot"] = NaN; rec["$(s)_its"] = NaN; rec["$(s)_inj"] = 0.0
                end
            end
            push!(rows, rec)
        end
    end
    return rows
end

# max injected level across the sparsify columns on one draw
row_inj(r) = maximum(Float64[get(r, "$(s)_inj", 0.0) for s in SPARSIFY]; init = 0.0)

# per-draw speedup base/target (base slower -> ratio > 1 -> target faster)
function speedups(rows, base, target)
    out = Float64[]
    for r in rows
        x = get(r, "$(base)_tot", NaN); y = get(r, "$(target)_tot", NaN)
        (isfinite(x) && isfinite(y) && y > 0) && push!(out, x / y)
    end
    return out
end

function report(rows, solvers)
    fams = unique(String[r["family"] for r in rows])
    sort!(fams)
    for fam in fams
        frows = [r for r in rows if r["family"] == fam]
        stalled = [r for r in frows if row_inj(r) > 0]
        clean = [r for r in frows if row_inj(r) == 0]
        println("\n=== $(fam)  ($(length(frows)) draws: $(length(stalled)) stalled, $(length(clean)) clean) ===")
        isempty(stalled) && (println("  no stalled draws (sparsify never injected)"); continue)

        injs = Float64[row_inj(r) for r in stalled]
        @printf("  injected levels on stalled draws: min %d, median %g, max %d\n",
                Int(minimum(injs)), median(injs), Int(maximum(injs)))

        println("  STALLED draws — median total s (median its):")
        for s in solvers
            t = medf(fin(Float64[r["$(s)_tot"] for r in stalled]))
            it = medf(fin(Float64[r["$(s)_its"] for r in stalled]))
            @printf("    %-16s %8s  (%s)\n", s, fmt(t), fmt(it))
        end
        println("  STALLED draws — per-draw median speedup (>1 = faster):")
        for s in SPARSIFY
            s in solvers || continue
            vk = speedups(stalled, REF, s)   # vs cmg-k-elim
            va = "ac" in solvers ? speedups(stalled, "ac", s) : Float64[]
            @printf("    %-16s vs %s %6s   vs ac %6s\n", s, REF,
                    isempty(vk) ? "--" : fmt(median(vk)),
                    isempty(va) ? "--" : fmt(median(va)))
        end
        # worst stalled draw for the reference (where cmg-k-elim is slowest vs ac)
        worst = nothing; wr = -Inf
        for r in stalled
            a = get(r, "ac_tot", NaN); k = get(r, "$(REF)_tot", NaN)
            (isfinite(a) && isfinite(k) && a > 0) || continue
            ratio = k / a
            ratio > wr && (wr = ratio; worst = r)
        end
        if worst !== nothing
            @printf("  worst stalled draw for %s (%s/ac = %s):\n", REF, REF, fmt(wr))
            @printf("    %s  n=%d\n", worst["matrix"], Int(worst["nv"]))
            for s in solvers
                @printf("      %-16s %8s s  (%s its, inj %d)\n", s,
                        fmt(get(worst, "$(s)_tot", NaN)), fmt(get(worst, "$(s)_its", NaN)),
                        Int(get(worst, "$(s)_inj", 0.0)))
            end
        end

        if !isempty(clean)
            println("  CLEAN draws — median total s (for contrast):")
            for s in solvers
                t = medf(fin(Float64[r["$(s)_tot"] for r in clean]))
                @printf("    %-16s %8s s\n", s, fmt(t))
            end
        end
    end
end

function main()
    inputs = String[]; solvers = copy(DEFAULT_SOLVERS)
    let i = 1
        while i <= length(ARGS)
            if ARGS[i] == "--solvers" && i < length(ARGS)
                solvers = String.(split(ARGS[i+1], ',')); i += 2
            else
                push!(inputs, ARGS[i]); i += 1
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
            @warn "no such path" path = inp
        end
    end
    isempty(files) && (println("no .jld2 files found"); return)
    rows = load_rows(files, solvers)
    println("loaded $(length(rows)) per-sample draws from $(length(files)) files")
    println("solvers: ", join(solvers, ", "))
    report(rows, solvers)
end

main()
