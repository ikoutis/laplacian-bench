#==========================================================
Prefetch the downloadable benchmark matrices into matrix-files/.

    julia --project=.. download_data.jl [--scale smoke|medium|paper]

Run this once on a machine WITH internet access (e.g. the Wulver login node)
before submitting jobs: compute nodes are assumed offline, and the runner's
loaders only fall back to downloading when a file is missing.

Fetches (idempotent — cached files are kept and skipped):
  * the curated SuiteSparse selection (matrix-files/suitesparse-selected.jld2),
  * the chimeraIPM matrices  uc.i<i>.eps<eps>.<cnt>  (SuiteSparse),
  * the spielmanIPM matrices sk<k>i<i>               (SuiteSparse),
and reports which SPE .mm files are present. The SPE matrices are not
downloadable from SuiteSparse — fetch them manually as described in
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

# ---- chimeraIPM -----------------------------------------------------------
is, js, cnts = scale === :paper  ? (1:5, 1:6, 1:6) :
               scale === :medium ? (1:2, 1:2, 1:3) :
                                   (1:1, 1:1, 1:1)
println("== chimeraIPM (uc.i*.eps*.*) ==")
for i in is, j in js
    curTargetEps = 1 / 10^j
    for cnt in cnts
        name = "uc.i$(i).eps$(curTargetEps).$(cnt)"
        if ipmAvailable(name)
            println("  $(name) already cached")
            global ok += 1
            continue
        end
        print("  $(name) ... ")
        M = loadSSCached(name)
        if M === nothing
            # matches downloadIPM.jl: not every (eps, cnt) slot exists upstream
            println("not available upstream; stopping this eps series")
            break
        end
        println("ok")
        global ok += 1
        M = nothing
        GC.gc()
    end
end

# ---- spielmanIPM ----------------------------------------------------------
ks, sis = scale === :paper  ? (100:100:600, 1:11) :
          scale === :medium ? (100:100:200, 1:3) :
                              (100:100:100, 1:1)
println("== spielmanIPM (sk<k>i<i>) ==")
for k in ks
    for i in sis
        name = "sk$(k)i$(i)"
        if ipmAvailable(name)
            println("  $(name) already cached")
            global ok += 1
            continue
        end
        print("  $(name) ... ")
        M = loadSSCached(name)
        if M === nothing
            println("not available upstream; stopping this k series")
            break
        end
        println("ok")
        global ok += 1
        M = nothing
        GC.gc()
    end
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
