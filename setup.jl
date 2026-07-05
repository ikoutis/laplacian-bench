#=
One-time environment setup for the chol-vs-kcycle benchmark.

    julia --project=. setup.jl

CombinatorialMultigrid's registered release (0.2.0) predates the K-cycle
solver (`cmg_solve`), so the package is added straight from GitHub, pinned to
a revision known to work with these scripts.

Overrides:
  CMG_DEV=/path/to/CombinatorialMultigrid.jl   use a local checkout (Pkg.develop)
  CMG_REV=<git rev>                            pin a different revision/branch
=#

using Pkg

Pkg.activate(@__DIR__)

const CMG_DEFAULT_REV = "08c515ed76aabc4caef08abc8bc60b98a549ea1e"

if haskey(ENV, "CMG_DEV")
    @info "Using local CombinatorialMultigrid checkout" path = ENV["CMG_DEV"]
    Pkg.develop(path = ENV["CMG_DEV"])
else
    rev = get(ENV, "CMG_REV", CMG_DEFAULT_REV)
    @info "Adding CombinatorialMultigrid from GitHub" rev
    Pkg.add(url = "https://github.com/ikoutis/CombinatorialMultigrid.jl", rev = rev)
end

Pkg.instantiate()
Pkg.precompile()

@info "Environment ready. Commit the generated Manifest.toml to pin it."
