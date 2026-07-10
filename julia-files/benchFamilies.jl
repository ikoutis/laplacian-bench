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
# The Dropbox ipmMat.zip unzips into this subdirectory; the IPM readers look
# here too so the sweep works whether it was flattened into matrix-files/ or
# left inside ipmMat/.
const IPM_SUBDIR = joinpath(MATRIX_DIR, "ipmMat")
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

"Read a MatrixMarket file from matrix-files/ (or matrix-files/ipmMat/); nothing if absent."
function loadMMCached(name)
    path = joinpath(MATRIX_DIR, name * ".mm")
    if !isfile(path)
        alt = joinpath(IPM_SUBDIR, name * ".mm")
        isfile(alt) || return nothing
        path = alt
    end
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
IPM matrices (chimeraIPM / spielmanIPM) are hosted in the SuiteSparse FlowIPM22
group. Prefer the .mm files the original scripts read, then a locally cached
.mat, then a live SuiteSparse download. Local files are cached/looked up under
the bare name (loadSSCached/ensureSSCached key the cache on the basename), so
only the download URL carries the FlowIPM22/ group prefix.
=#
ipmSSName(name) = "FlowIPM22/" * name

# FlowIPM22 stores each graph as "Undirected Weighted Graph" (adjacency, zero
# diagonal); older IPM files were Laplacians (positive diagonal). Return the
# adjacency for either: a Laplacian goes through adj(), an adjacency is used
# as-is. This keeps the :lap path (which rebuilds lap(adjacency)) correct and
# avoids handing CMG a positive-off-diagonal matrix.
ipmAdjacency(M) = all(iszero, diag(M)) ? M : adj(M)[1]

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
    return loadSSCached(ipmSSName(name); download = download)
end

ipmAvailable(name) = isfile(joinpath(MATRIX_DIR, name * ".mm")) ||
                     isfile(joinpath(MATRIX_DIR, name * ".mat"))

# True if an IPM matrix is (or can be made) available: already local, or
# downloadable now (and downloads allowed). Used for instance-list discovery
# — like the *_ac.jl scripts, the sweep stops at the first name that isn't
# obtainable. Downloads eagerly so the run-time loader reads from cache.
ensureIPMCached(name; download = ALLOW_DOWNLOAD[]) =
    ipmAvailable(name) || (ensureSSCached(ipmSSName(name); download = download) !== nothing)

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

# Natural-order sort key: compare the sequence of embedded integers first (so
# sk100i2 < sk100i10), then the raw string.
natkey(s::AbstractString) = ([parse(Int, m.match) for m in eachmatch(r"\d+", s)], String(s))

# Bare names (no ".mm") of locally-staged IPM Matrix Market files whose bare
# name matches `rx`, naturally sorted. These come from the manual ipmMat.zip
# (the paper's full IPM sweep), which is NOT auto-downloadable. Scans both
# matrix-files/ and matrix-files/ipmMat/ (the zip's own subdirectory), so the
# sweep is found whether it was flattened or left inside ipmMat/.
function localIPMNames(rx::Regex)
    names = String[]
    for dir in (MATRIX_DIR, IPM_SUBDIR)
        isdir(dir) || continue
        for f in readdir(dir)
            endswith(f, ".mm") || continue
            bare = f[1:end-3]
            occursin(rx, bare) && push!(names, bare)
        end
    end
    unique!(names)          # a name in both dirs -> once (top level wins on load)
    sort!(names; by = natkey)
    return names
end

ipmInst(fam, nm) = BenchInstance("$(fam) $(nm)",
    (baseseed, rep) -> begin
        L = loadMMCached(nm)
        L === nothing && return nothing
        (:lap, ipmAdjacency(L), nm)
    end)

# chimeraIPM. If the manual full sweep is staged (in matrix-files/ or its
# ipmMat/ subdir), use the paper's full set (capped by scale on the chimera
# index i); otherwise fall back to the auto-downloadable FlowIPM22 subset
# uni_chimera_i<i>. The Dropbox ipmMat.zip names them
# uni_chimera.n<n>.i<i>.eps<eps>.<cnt>.mm; the short uc.i<i>.eps<eps>.<cnt> form
# is also accepted.
function chimeraIPMInstances(scale; n = nothing, limit = nothing)
    imax = scale === :paper ? 5 : scale === :medium ? 3 : 1
    full = localIPMNames(r"^(?:uni_chimera\.n\d+|uc)\.i\d+\.eps")
    if !isempty(full)
        keep = filter(nm -> (m = match(r"\.i(\d+)\.eps", nm);
                             m !== nothing && parse(Int, m.captures[1]) <= imax), full)
        insts = [ipmInst("chimeraIPM", nm) for nm in keep]
        isempty(insts) && @warn "chimeraIPM: staged chimera-IPM .mm present but none with i <= $(imax)"
        return applyLimit(insts, limit)
    end
    insts = BenchInstance[]
    for i in 1:imax
        name = "uni_chimera_i$(i)"
        ensureIPMCached(name) || continue
        push!(insts, BenchInstance("chimeraIPM $(name)",
            (baseseed, rep) -> begin
                L = loadIPM(name)
                L === nothing && return nothing
                (:lap, ipmAdjacency(L), name)
            end))
    end
    if isempty(insts)
        @warn "chimeraIPM: no chimera-IPM .mm (manual ipmMat.zip: uni_chimera.n*.i*.eps*.mm) and no uni_chimera_i*.mm/.mat — run download_data.jl or stage ipmMat.zip"
    end
    return applyLimit(insts, limit)
end

# spielmanIPM. If the manual full sweep is staged (in matrix-files/ or its
# ipmMat/ subdir), use the paper's full set (k = 100..600 by scale, i = 1..11);
# otherwise fall back to the auto FlowIPM22 subset Spielman_k<k>. The Dropbox
# ipmMat.zip names them spielman.k<k>.low<lo>.up<up>.i<i>.mm; the short sk<k>i<i>
# form is also accepted. NOTE the largest instances are big-memory (this family
# runs in the big-memory tier; see BIG_FIXED in run_chol_vs_kcycle.sh); cap with
# --limit to skip the biggest.
function spielmanIPMInstances(scale; n = nothing, limit = nothing)
    kmax = scale === :paper ? 600 : scale === :medium ? 300 : 100
    full = localIPMNames(r"^(?:spielman\.k\d+\..*\.i\d+|sk\d+i\d+)$")
    if !isempty(full)
        keep = filter(nm -> (m = match(r"^(?:spielman\.k|sk)(\d+)", nm);
                             m !== nothing && parse(Int, m.captures[1]) <= kmax), full)
        insts = [ipmInst("spielmanIPM", nm) for nm in keep]
        isempty(insts) && @warn "spielmanIPM: staged spielman-IPM .mm present but none with k <= $(kmax)"
        return applyLimit(insts, limit)
    end
    insts = BenchInstance[]
    for k in 100:100:kmax
        name = "Spielman_k$(k)"
        ensureIPMCached(name) || continue
        push!(insts, BenchInstance("spielmanIPM $(name)",
            (baseseed, rep) -> begin
                L = loadIPM(name)
                L === nothing && return nothing
                (:lap, ipmAdjacency(L), name)
            end))
    end
    if isempty(insts)
        @warn "spielmanIPM: no spielman-IPM .mm (manual ipmMat.zip: spielman.k*.low*.up*.i*.mm) and no Spielman_k*.mm/.mat — run download_data.jl or stage ipmMat.zip"
    end
    return applyLimit(insts, limit)
end

# Chimera families: the paper reports statistics over C instances per vertex
# count, the Chimeras with seed indices i = 1..C (chimeraAndIPM.tex; "we always
# choose the first seed indices ... we do not exclude any Chimeras"). C shrinks
# with size because the large graphs are expensive. Each (size, i) is one
# BenchInstance; `reps` then repeats the *timing* of these same instances rather
# than drawing new graphs (unlike the earlier one-per-size design).
const CHIMERA_COUNTS_PAPER = Dict(10^4 => 103, 10^5 => 105, 10^6 => 23, 10^7 => 8)
chimeraCount(scale, sz) =
    scale === :paper  ? get(CHIMERA_COUNTS_PAPER, Int(sz), 1) :
    scale === :medium ? 5 :
                        2

# Helper so each closure captures its own (sz, i) — a bare for-loop closure over
# the loop variables is a classic Julia capture footgun.
chimeraInstance(gen, genname::String, kind::Symbol, sz, i) =
    BenchInstance("$(genname) n=$(sz) i=$(i)",
        (baseseed, rep) -> (kind, gen(Int64(sz), i), "$(genname)($(Int64(sz)),$(i))"))

function chimeraInstances(gen, genname::String, kind::Symbol, scale; n = nothing, limit = nothing)
    sizes = n !== nothing ? n :
            scale === :paper  ? [10^4, 10^5, 10^6, 10^7] :
            scale === :medium ? [10^4, 10^5] :
                                [10^4]
    insts = BenchInstance[]
    for sz in sizes, i in 1:chimeraCount(scale, sz)
        push!(insts, chimeraInstance(gen, genname, kind, sz, i))
    end
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
