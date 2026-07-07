#=
One-time environment setup for the chol-vs-kcycle benchmark.

    julia --project=. setup.jl

CombinatorialMultigrid's registered release (0.2.0) predates the K-cycle
solver (`cmg_solve`), so the package is added straight from GitHub, pinned to
a revision known to work with these scripts.

Overrides:
  CMG_DEV=/path/to/CombinatorialMultigrid.jl   use a local checkout (Pkg.develop)
  CMG_REV=<git rev>                            pin a different revision/branch
  SETUP_SKIP_PRECOMPILE=1                       resolve/download only, no precompile
                                               (precompilation then happens on
                                               first `using` — useful when this
                                               step runs where precompiling is
                                               restricted, e.g. a cluster login
                                               node)

Resolution order for CombinatorialMultigrid: CMG_DEV, then an in-repo checkout
at ./CombinatorialMultigrid.jl (git-ignored; see performance-experiments/wulver/
fetch_cmg.sh), then CMG_REV / the default GitHub revision. The in-repo checkout
lets the whole environment come from one `git clone` with no run-time GitHub
fetch — handy on clusters whose compute nodes have no internet.
=#

using Pkg

Pkg.activate(@__DIR__)

const CMG_DEFAULT_REV = "08c515ed76aabc4caef08abc8bc60b98a549ea1e"
const VENDORED_CMG = joinpath(@__DIR__, "CombinatorialMultigrid.jl")

if haskey(ENV, "CMG_DEV")
    @info "Using local CombinatorialMultigrid checkout" path = ENV["CMG_DEV"]
    Pkg.develop(path = ENV["CMG_DEV"])
elseif !haskey(ENV, "CMG_REV") && isdir(VENDORED_CMG)
    @info "Using in-repo CombinatorialMultigrid checkout" path = VENDORED_CMG
    Pkg.develop(path = VENDORED_CMG)
else
    rev = get(ENV, "CMG_REV", CMG_DEFAULT_REV)
    @info "Adding CombinatorialMultigrid from GitHub" rev
    Pkg.add(url = "https://github.com/ikoutis/CombinatorialMultigrid.jl", rev = rev)
end

Pkg.instantiate()
if get(ENV, "SETUP_SKIP_PRECOMPILE", "0") == "1"
    @info "Skipping Pkg.precompile() (SETUP_SKIP_PRECOMPILE=1); it runs on first `using`."
else
    Pkg.precompile()
end

@info "Environment ready. Commit the generated Manifest.toml to pin it."
