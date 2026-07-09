#==========================================================
Prefetch the downloadable benchmark matrices into matrix-files/.

    julia --project=.. download_data.jl [--scale smoke|medium|paper]

Run this once on a machine WITH internet access (e.g. the Wulver login node)
before submitting jobs: compute nodes are assumed offline, and the runner's
loaders only fall back to downloading when a file is missing.

Fetches (idempotent — cached files are kept and skipped):
  * the curated SuiteSparse selection (matrix-files/suitesparse-selected.jld2),
  * the chimeraIPM  FlowIPM22 SUBSET  uni_chimera_i<i>  (SuiteSparse),
  * the spielmanIPM FlowIPM22 SUBSET  Spielman_k<k>     (SuiteSparse),
and reports which SPE .mm files are present.

Two data sources are MANUAL (not fetched here):
  * The paper's FULL IPM sweep — uc.i<i>.eps<eps>.<cnt>.mm (chimera-IPM) and
    sk<k>i<i>.mm (Spielman-IPM) — ships as ipmMat.zip (Tutorial.md). Unzip it
    into matrix-files/ and the loaders prefer it over the FlowIPM22 subset
    fetched here (see benchFamilies.jl chimeraIPMInstances/spielmanIPMInstances).
  * SPE is not downloadable from SuiteSparse — fetch it manually as described in
    Tutorial.md ("SPE benchmark") and drop spe{0.5m,2m,4m,8m,16m}.mm into
    matrix-files/.
===========================================================#

using JLD2

include(joinpath(@__DIR__, "..", "julia-files", "benchFamilies.jl"))

scale = :paper
let args = ARGS
    i = 1
    while i <= length(args)
        if args[i] == "--scale" && i < length(args)
            global scale = Symbol(args[i+1])
            i += 1
        else
            println(stderr, "usage: julia --project=.. download_data.jl [--scale smoke|medium|paper]")
            exit(1)
        end
        i += 1
    end
end
scale in (:smoke, :medium, :paper) || (println(stderr, "bad --scale $(scale)"); exit(1))

allow_download!(true)   # this is the prefetch tool; always permit downloads

println("prefetching benchmark data at scale = $(scale) into $(MATRIX_DIR)")
ok, failed = 0, 0

# ---- SuiteSparse selection ------------------------------------------------
names = suitesparseNames()
if scale === :medium
    names = names[1:min(length(names), 6)]
elseif scale === :smoke
    names = names[1:min(length(names), 1)]
end
println("== suitesparse selection: $(length(names)) matrices ==")
for name in names
    print("  $(name) ... ")
    M = loadSSCached(name)
    if M === nothing
        println("FAILED")
        global failed += 1
    else
        println("ok ($(size(M,1)) x $(size(M,2)))")
        global ok += 1
    end
    M = nothing
    GC.gc()
end

# ---- chimeraIPM (FlowIPM22 uni_chimera_i*) --------------------------------
imax = scale === :paper ? 5 : scale === :medium ? 3 : 1
println("== chimeraIPM (FlowIPM22/uni_chimera_i*) ==")
for i in 1:imax
    name = "uni_chimera_i$(i)"
    if ipmAvailable(name)
        println("  $(name) already cached")
        global ok += 1
        continue
    end
    print("  $(name) ... ")
    if ensureIPMCached(name)
        println("ok")
        global ok += 1
    else
        println("not available upstream")
        global failed += 1
    end
    GC.gc()
end

# ---- spielmanIPM (FlowIPM22 Spielman_k*; k500/k600 are multi-GB) ----------
ks = scale === :paper  ? (100:100:600) :
     scale === :medium ? (100:100:300) :
                         (100:100:100)
println("== spielmanIPM (FlowIPM22/Spielman_k*) ==")
for k in ks
    name = "Spielman_k$(k)"
    if ipmAvailable(name)
        println("  $(name) already cached")
        global ok += 1
        continue
    end
    print("  $(name) ... ")
    if ensureIPMCached(name)
        println("ok")
        global ok += 1
    else
        println("not available upstream")
        global failed += 1
    end
    GC.gc()
end

# ---- SPE (manual) -----------------------------------------------------------
spe_sizes = scale === :paper  ? ["0.5m", "2m", "4m", "8m", "16m"] :
            scale === :medium ? ["0.5m", "2m"] :
                                ["0.5m"]
println("== spe (manual download) ==")
missing_spe = String[]
for sz in spe_sizes
    if isfile(joinpath(MATRIX_DIR, "spe$(sz).mm"))
        println("  spe$(sz).mm present")
    else
        println("  spe$(sz).mm MISSING")
        push!(missing_spe, sz)
    end
end
if !isempty(missing_spe)
    println("  -> the spe family will skip these sizes; see Tutorial.md for the SPE download link")
end

println("done: $(ok) fetched/cached, $(failed) failures")
exit(failed == 0 ? 0 : 1)
