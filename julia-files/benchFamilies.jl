#==========================================================
Registry of the benchmark problem families for the chol-vs-kcycle runner.

Each family reproduces the instance sweep of its `*_ac.jl` experiment script
(the "paper" scale) and adds smaller "medium"/"smoke" sweeps so the whole
suite can be exercised quickly. Generator calls are lifted verbatim from:

  aniso_ac.jl, checkered_ac.jl, uniform_grid_ac.jl, wgrid_ac.jl,
  sachdeva_star_ac.jl, suitesparse_ac.jl, spe_ac.jl,
  chimeraIPM_ac.jl, spielmanIPM_ac.jl,
  uni_chimera_ac.jl, uni_bndry_chimera_ac.jl,
  wted_chimera_ac.jl, wted_bndry_chimera_ac.jl

An instance's `load(baseseed, rep)` returns `(kind, mat, label)` where
`kind` is `:lap` (mat = adjacency matrix) or `:sddm` (mat = SDDM matrix),
or `nothing` when the instance's data is unavailable (missing download);
the runner skips those with a warning and continues.

Chimera instances consume `(baseseed, rep)` as the chimera index
`i = 1000*baseseed + rep`, so different seeds give disjoint graph streams
and reruns with the same (seed, rep) are identical. Deterministic families
ignore the arguments (their repetition variance comes from the seeded RHS
and the solvers' internal randomness).

All file access is anchored at this file's location, never the cwd, so the
runner works from any directory (as Slurm requires). Downloads cache into
matrix-files/ and are never deleted.
===========================================================#

using Laplacians
using SparseArrays
using LinearAlgebra
using JLD2
using MAT
using MatrixMarket
using Downloads

include(joinpath(@__DIR__, "isLap.jl"))

# star_join comes from Laplacians (as sachdeva_star_ac.jl assumes); fall back
# to the local copy for Laplacians versions that lack it.
if !isdefined(Laplacians, :star_join)
    include(joinpath(@__DIR__, "star_join.jl"))
end

const MATRIX_DIR = normpath(joinpath(@__DIR__, "..", "matrix-files"))
const SS_URL = "http://sparse-files.engr.tamu.edu/mat/"

# Global download policy, set once by the runner (chol_vs_kcycle.jl) or
# download_data.jl. When true (default), every downloadable family fetches
# missing matrices on demand and caches them under matrix-files/; when false
# (runner --offline / CVK_OFFLINE), no run-time family touches the network and
# only already-cached/prefetched data is used. Applied uniformly to the
# SuiteSparse and IPM families (SPE is Dropbox-only and always local).
const ALLOW_DOWNLOAD = Ref(true)
allow_download!(b::Bool) = (ALLOW_DOWNLOAD[] = b)

struct BenchInstance
    base::String        # short name for logs/progress
    load::Function      # (baseseed, rep) -> nothing | (kind, mat, label)
end

struct BenchFamily
    name::String
    tol::Float64        # default tolerance (spe uses 5e-9 like spe_ac.jl)
    sized::Bool         # true: accepts a --n size list (chimera families)
    instances::Function # (scale::Symbol; n=nothing, limit=nothing) -> Vector{BenchInstance}
end

# ------------------------------------------------------------------ loaders

"Read a MatrixMarket file from matrix-files/; nothing if absent."
function loadMMCached(name)
    path = joinpath(MATRIX_DIR, name * ".mm")
    isfile(path) || return nothing
    return MatrixMarket.mmread(path)
end

"Read Problem.A from a SuiteSparse-format .mat file."
function readSSMat(path)
    file = MAT.matopen(path)
    M = try
        read(file, "Problem")["A"]
    finally
        close(file)
    end
    return M
end

#=
Ensure a SuiteSparse .mat is cached under matrix-files/ (keyed by the base
name, as downloadSS names it) and never deleted. Downloads it if missing and
downloads are allowed. Returns the local path, or nothing if unavailable
(offline and not cached, or the fetch failed). Does not parse the matrix.
=#
function ensureSSCached(fullname; download = ALLOW_DOWNLOAD[])
    base = String(split(fullname, '/')[end])
    path = joinpath(MATRIX_DIR, base * ".mat")
    isfile(path) && return path
    download || return nothing
    try
        tmp = path * ".part"
        Downloads.download(string(SS_URL, fullname, ".mat"), tmp)
        mv(tmp, path; force = true)
        return path
    catch e
        @warn "download failed for SuiteSparse matrix $(fullname)" exception = e
        return nothing
    end
end

"Cached SuiteSparse fetch: returns the matrix, or nothing if unavailable."
function loadSSCached(fullname; download = ALLOW_DOWNLOAD[])
    path = ensureSSCached(fullname; download = download)
    path === nothing && return nothing
    try
        return readSSMat(path)
    catch e
        @warn "could not read $(path)" exception = e
        return nothing
    end
end

#=
IPM matrices (chimeraIPM / spielmanIPM): prefer the .mm files the original
scripts read, then a locally cached .mat (as downloadIPM.jl fetches them),
then a live SuiteSparse download.
=#
function loadIPM(name; download = ALLOW_DOWNLOAD[])
    M = loadMMCached(name)
    M === nothing || return M
    path = joinpath(MATRIX_DIR, name * ".mat")
    if isfile(path)
        try
            return readSSMat(path)
        catch e
            @warn "could not read $(path)" exception = e
        end
    end
    return loadSSCached(name; download = download)
end

ipmAvailable(name) = isfile(joinpath(MATRIX_DIR, name * ".mm")) ||
                     isfile(joinpath(MATRIX_DIR, name * ".mat"))

# True if an IPM matrix is (or can be made) available: already local, or
# downloadable now (and downloads allowed). Used for instance-list discovery
# — like the *_ac.jl scripts, the sweep stops at the first name that isn't
# obtainable. Downloads eagerly so the run-time loader reads from cache.
ensureIPMCached(name; download = ALLOW_DOWNLOAD[]) =
    ipmAvailable(name) || (ensureSSCached(name; download = download) !== nothing)

"suitesparse-selected.jld2 holds the curated SuiteSparse matrix name list."
function suitesparseNames()
    path = joinpath(MATRIX_DIR, "suitesparse-selected.jld2")
    if !isfile(path)
        @warn "missing $(path); suitesparse family is empty"
        return String[]
    end
    try
        return JLD2.load(path, "names")
    catch e
        @warn "could not read $(path); suitesparse family is empty" exception = e
        return String[]
    end
end

applyLimit(v, limit) = limit === nothing ? v : v[1:min(length(v), Int(limit))]

# --------------------------------------------------------- family builders

# Deterministic SDDM generator families (aniso/checkered/uniform_grid/wgrid).
sddmInstance(base, label, gen) =
    BenchInstance(base, (baseseed, rep) -> (:sddm, gen(), label))

function anisoInstances(scale; n = nothing, limit = nothing)
    nnz_var, xis = scale === :paper  ? (2e8, [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]) :
                   scale === :medium ? (2e6, [1e-2, 1, 1e2]) :
                                       (2e5, [1e-2, 1e2])
    insts = [sddmInstance("aniso xi=$(xi)", "aniso_grid_sddm($(nnz_var),$(xi))",
                          () -> aniso_grid_sddm(Int64(nnz_var), xi))
             for xi in xis]
    return applyLimit(insts, limit)
end

function checkeredInstances(scale; n = nothing, limit = nothing)
    w = 1e7
    nnz_var, blocks = scale === :paper  ? (2e8, [2, 4, 8, 16, 32, 64, 128]) :
                      scale === :medium ? (2e6, [2, 8]) :
                                          (2e5, [2])
    insts = [sddmInstance("checkered b=$(blk)",
                          "checkered_grid_sddm($(nnz_var), $(blk), $(blk), $(blk), $(w))",
                          () -> checkered_grid_sddm(Int64(nnz_var), blk, blk, blk, w))
             for blk in blocks]
    return applyLimit(insts, limit)
end

function uniformGridInstances(scale; n = nothing, limit = nothing)
    nnzs = scale === :paper  ? [2e6, 2e7, 2e8] :
           scale === :medium ? [2e6] :
                               [2e5]
    insts = [sddmInstance("uniform_grid nnz=$(nnz_var)", "uniform_grid_sddm($(nnz_var))",
                          () -> uniform_grid_sddm(Int64(nnz_var)))
             for nnz_var in nnzs]
    return applyLimit(insts, limit)
end

function wgridInstances(scale; n = nothing, limit = nothing)
    nnz_var, ws = scale === :paper  ? (2e8, [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]) :
                  scale === :medium ? (2e6, [1e-2, 1, 1e2]) :
                                      (2e5, [1e-2, 1e2])
    insts = [sddmInstance("wgrid w=$(w)", "wgrid_sddm($(nnz_var), $(w))",
                          () -> wgrid_sddm(Int64(nnz_var), w))
             for w in ws]
    return applyLimit(insts, limit)
end

function sachdevaStarInstances(scale; n = nothing, limit = nothing)
    ks = scale === :paper  ? collect(100:50:800) :
         scale === :medium ? collect(100:100:300) :
                             [100]
    insts = [BenchInstance("sachdeva_star k=$(k)",
                 (baseseed, rep) -> begin
                     l = Int64(k / 2)
                     A = star_join(complete_graph(k), l)
                     (:lap, A, "star_join(complete_graph($(k)), $(l))")
                 end)
             for k in ks]
    return applyLimit(insts, limit)
end

function suitesparseInstances(scale; n = nothing, limit = nothing)
    names = suitesparseNames()
    if scale === :medium
        names = names[1:min(length(names), 6)]
    elseif scale === :smoke
        names = names[1:min(length(names), 1)]
    end
    insts = [BenchInstance("suitesparse $(name)",
                 (baseseed, rep) -> begin
                     M = loadSSCached(name)
                     M === nothing && return nothing
                     if isLap(M)
                         a, _ = adj(M)
                         return (:lap, a, String(name))
                     else
                         return (:sddm, M, String(name))
                     end
                 end)
             for name in names]
    return applyLimit(insts, limit)
end

function speInstances(scale; n = nothing, limit = nothing)
    sizes = scale === :paper  ? ["0.5m", "2m", "4m", "8m", "16m"] :
            scale === :medium ? ["0.5m", "2m"] :
                                ["0.5m"]
    insts = [BenchInstance("spe$(sz)",
                 (baseseed, rep) -> begin
                     mat = loadMMCached("spe$(sz)")
                     if mat === nothing
                         @warn "spe$(sz).mm not found in matrix-files/ — skipping (SPE files must be fetched manually; see Tutorial.md)"
                         return nothing
                     end
                     (:sddm, -mat, "spe$(sz)")   # stored negative-definite, as in spe_ac.jl
                 end)
             for sz in sizes]
    return applyLimit(insts, limit)
end

# chimeraIPM: SuiteSparse-hosted Laplacians named uc.i<i>.eps<eps>.<cnt>.
# The original script iterates cnt until the first missing file; we enumerate
# availability the same way (checking .mm and cached .mat).
function chimeraIPMInstances(scale; n = nothing, limit = nothing)
    is, js, cnts = scale === :paper  ? (1:5, 1:6, 1:6) :
                   scale === :medium ? (1:2, 1:2, 1:3) :
                                       (1:1, 1:1, 1:1)
    insts = BenchInstance[]
    for i in is, j in js
        curTargetEps = 1 / 10^j
        for cnt in cnts
            name = "uc.i$(i).eps$(curTargetEps).$(cnt)"
            if !ensureIPMCached(name)
                break   # matches chimeraIPM_ac.jl: stop at the first missing count
            end
            push!(insts, BenchInstance("chimeraIPM $(name)",
                (baseseed, rep) -> begin
                    L = loadIPM(name)
                    L === nothing && return nothing
                    a, _ = adj(L)
                    # the original reused a stale testName here; label properly
                    (:lap, a, name)
                end))
        end
    end
    if isempty(insts)
        @warn "chimeraIPM: no uc.i*.eps*.* files under matrix-files/ — run performance-experiments/download_data.jl first"
    end
    return applyLimit(insts, limit)
end

function spielmanIPMInstances(scale; n = nothing, limit = nothing)
    ks, is = scale === :paper  ? (100:100:600, 1:11) :
             scale === :medium ? (100:100:200, 1:3) :
                                 (100:100:100, 1:1)
    insts = BenchInstance[]
    for k in ks, i in is
        name = "sk$(k)i$(i)"
        if !ensureIPMCached(name)
            break   # matches spielmanIPM_ac.jl: stop at the first missing index
        end
        push!(insts, BenchInstance("spielmanIPM $(name)",
            (baseseed, rep) -> begin
                L = loadIPM(name)
                L === nothing && return nothing
                a, _ = adj(L)
                (:lap, a, name)
            end))
    end
    if isempty(insts)
        @warn "spielmanIPM: no sk*i*.mm/.mat files under matrix-files/ — run performance-experiments/download_data.jl first"
    end
    return applyLimit(insts, limit)
end

# Chimera families: one instance per size; the chimera index is derived from
# (baseseed, rep) so instances are reproducible and seed streams are disjoint.
chimeraIndex(baseseed, rep) = 1000 * Int(baseseed) + Int(rep)

function chimeraInstances(gen, genname::String, kind::Symbol, scale; n = nothing, limit = nothing)
    sizes = n !== nothing ? n :
            scale === :paper  ? [10^4, 10^5, 10^6, 10^7] :
            scale === :medium ? [10^4, 10^5] :
                                [10^4]
    insts = [BenchInstance("$(genname) n=$(sz)",
                 (baseseed, rep) -> begin
                     i = chimeraIndex(baseseed, rep)
                     mat = gen(Int64(sz), i)
                     (kind, mat, "$(genname)($(Int64(sz)),$(i))")
                 end)
             for sz in sizes]
    return applyLimit(insts, limit)
end

# Laplacians 1.4's `uni_bndry_chimera` / `wted_bndry_chimera` build the interior
# index set with a Float64 step range (`setdiff(1:n, 1:n^(1/3):n)`), so `int`
# becomes a `Vector{Float64}` and `L[int, int]` throws `invalid index: 2.0 of
# type Float64` on modern Julia. Reimplement with the integer step Laplacians
# master uses (`ceil(Int, n^(1/3))`) so the boundary nodes are actually removed.
function uni_bndry_chimera_fixed(n::Integer, k::Integer)
    a = Laplacians.chimera(n, k)
    Laplacians.unweight!(a)
    L = Laplacians.lap(a)
    int = setdiff(1:n, 1:ceil(Int, n^(1 / 3)):n)
    return L[int, int]
end

function wted_bndry_chimera_fixed(n::Integer, k::Integer)
    a = Laplacians.wted_chimera(n, k)
    L = Laplacians.lap(a)
    int = setdiff(1:n, 1:ceil(Int, n^(1 / 3)):n)
    return L[int, int]
end

# ---------------------------------------------------------------- registry

const FAMILIES = Dict{String,BenchFamily}(
    "uniform_grid" => BenchFamily("uniform_grid", 1e-8, false, uniformGridInstances),
    "aniso"        => BenchFamily("aniso", 1e-8, false, anisoInstances),
    "wgrid"        => BenchFamily("wgrid", 1e-8, false, wgridInstances),
    "checkered"    => BenchFamily("checkered", 1e-8, false, checkeredInstances),
    "sachdeva_star" => BenchFamily("sachdeva_star", 1e-8, false, sachdevaStarInstances),
    "suitesparse"  => BenchFamily("suitesparse", 1e-8, false, suitesparseInstances),
    "spe"          => BenchFamily("spe", 5e-9, false, speInstances),
    "chimeraIPM"   => BenchFamily("chimeraIPM", 1e-8, false, chimeraIPMInstances),
    "spielmanIPM"  => BenchFamily("spielmanIPM", 1e-8, false, spielmanIPMInstances),
    "uni_chimera" => BenchFamily("uni_chimera", 1e-8, true,
        (scale; n = nothing, limit = nothing) ->
            chimeraInstances(uni_chimera, "uni_chimera", :lap, scale; n = n, limit = limit)),
    "uni_bndry_chimera" => BenchFamily("uni_bndry_chimera", 1e-8, true,
        (scale; n = nothing, limit = nothing) ->
            chimeraInstances(uni_bndry_chimera_fixed, "uni_bndry_chimera", :sddm, scale; n = n, limit = limit)),
    "wted_chimera" => BenchFamily("wted_chimera", 1e-8, true,
        (scale; n = nothing, limit = nothing) ->
            chimeraInstances(wted_chimera, "wted_chimera", :lap, scale; n = n, limit = limit)),
    "wted_bndry_chimera" => BenchFamily("wted_bndry_chimera", 1e-8, true,
        (scale; n = nothing, limit = nothing) ->
            chimeraInstances(wted_bndry_chimera_fixed, "wted_bndry_chimera", :sddm, scale; n = n, limit = limit)),
)

# Order used by `all`, the superscript and the Slurm manifests.
const FAMILY_ORDER = [
    "uniform_grid", "aniso", "wgrid", "checkered",
    "sachdeva_star", "suitesparse", "spe",
    "chimeraIPM", "spielmanIPM",
    "uni_chimera", "uni_bndry_chimera", "wted_chimera", "wted_bndry_chimera",
]
