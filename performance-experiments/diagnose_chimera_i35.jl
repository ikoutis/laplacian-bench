#!/usr/bin/env julia
# Dissect the chimera draw that defeats CMG: uni_chimera(100000, 35).
#
# In the paper-count run this draw (and its wted/bndry overlays) hit the 1000-
# iteration cap for every CMG variant while ApproxChol converged. This rebuilds
# the graph and reports (a) its structure, (b) the CMG coarsening hierarchy, and
# (c) a convergence sketch — relres at increasing iteration caps — for cmg-v,
# cmg-k, and cmg-k-elim, so we can see plateau vs slow-decrease. A known-good
# draw (i=1) at the same size is shown alongside for contrast.
#
# Usage (compute node, after `source wulver/env_wulver.sh`):
#   julia --project=.. diagnose_chimera_i35.jl [N] [BAD_I] [GOOD_I]
# Defaults: N=100000, BAD_I=35, GOOD_I=1.

using Laplacians, CombinatorialMultigrid, SparseArrays, LinearAlgebra, Random, Printf

N     = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 100_000
BAD_I = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 35
GOOD_I= length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1

# --- component sizes via iterative BFS on the adjacency ---
function component_sizes(a)
    n = size(a, 1); seen = falses(n); rows = rowvals(a); sizes = Int[]
    for s in 1:n
        seen[s] && continue
        cnt = 0; stack = [s]; seen[s] = true
        while !isempty(stack)
            u = pop!(stack); cnt += 1
            for idx in nzrange(a, u)
                v = rows[idx]
                if !seen[v]; seen[v] = true; push!(stack, v); end
            end
        end
        push!(sizes, cnt)
    end
    return sort(sizes; rev = true)
end

function structure(name, a)
    n = size(a, 1); deg = diff(a.colptr); w = nonzeros(a)
    comps = component_sizes(a)
    @printf("  %-22s n=%d  m=%d  avg_deg=%.2f\n", name, n, nnz(a) ÷ 2, sum(deg) / n)
    @printf("     degree:   min=%d  max=%d  deg1=%d  deg2=%d  (deg<=2: %.1f%%)\n",
            minimum(deg), maximum(deg), count(==(1), deg), count(==(2), deg),
            100 * count(<=(2), deg) / n)
    @printf("     weights:  min=%.3g  max=%.3g  ratio=%.3g\n",
            minimum(w), maximum(w), maximum(w) / minimum(w))
    @printf("     components: %d  (largest %s)\n", length(comps),
            join(comps[1:min(5, end)], ", "))
    # top-degree hubs
    p = partialsortperm(deg, 1:min(5, n); rev = true)
    @printf("     top degrees: %s\n", join(deg[p], ", "))
end

function hierarchy(name, a)
    local H
    try
        (_, H) = cmg_preconditioner_lap(lap(a); cycle = :kcycle, eliminate = false)
    catch e
        @printf("  %-22s BUILD ERROR: %s\n", name, sprint(showerror, e))
        return
    end
    ns = [Int(h.n) for h in H]
    ratios = [@sprintf("%.2f", ns[i+1] / ns[i]) for i in 1:length(ns)-1]
    @printf("  %-22s %d levels  sizes=%s\n", name, length(ns), join(ns, "->"))
    @printf("     coarsening ratios: %s\n", join(ratios, ", "))
end

function sketch(name, a)
    Random.seed!(1234)
    b = lap(a) * randn(size(a, 1)); b ./= norm(b)
    variants = [("cmg-v", :vcycle, false), ("cmg-k", :kcycle, false),
                ("cmg-k-elim", :kcycle, true)]
    caps = [50, 100, 300, 1000, 3000, 10000]
    println("  $(name):")
    for (lbl, cyc, elim) in variants
        local H
        try
            (_, H) = cmg_preconditioner_lap(lap(a); cycle = cyc, eliminate = elim)
        catch e
            println("    $(rpad(lbl, 11)) BUILD ERROR: $(sprint(showerror, e))")
            continue
        end
        line = "    $(rpad(lbl, 11))"
        for cap in caps
            _, st = cmg_solve(H, b; tol = 1e-8, maxit = cap, cycle = cyc)
            mark = st.converged ? "*" : " "
            line *= @sprintf(" %6d:%.1e%s", cap, st.relres, mark)
        end
        println(line)
    end
end

println("Building graphs …")
abad  = uni_chimera(N, BAD_I)
agood = uni_chimera(N, GOOD_I)

println("\n== Structure ==")
structure("uni_chimera($N,$BAD_I) BAD", abad)
structure("uni_chimera($N,$GOOD_I) good", agood)

println("\n== CMG coarsening hierarchy (eliminate=false) ==")
hierarchy("uni_chimera($N,$BAD_I) BAD", abad)
hierarchy("uni_chimera($N,$GOOD_I) good", agood)

println("\n== Convergence sketch: relres at cap  (* = converged to 1e-8) ==")
sketch("uni_chimera($N,$BAD_I) BAD", abad)
sketch("uni_chimera($N,$GOOD_I) good", agood)
