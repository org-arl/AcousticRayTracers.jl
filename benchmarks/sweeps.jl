# ---------------------------------------------------------------------------
# Sweep grids
#
# One-at-a-time sensitivity around defaults, per model. The last entry of each
# list is taken as the "finest" self-convergence reference.

"RaySolver option grids for a scenario (ds scaled to water depth)."
function raysolver_sweeps(scn)
  depth = -minimum(scn.grid.zrange)  # grid depth extent as proxy for water depth
  Dict(
    :nbeams => [(nbeams=n,) for n ∈ [0, 250, 500, 1000, 2000, 4000]],
    :ds => [(ds=d,) for d ∈ [depth/2, depth/5, 0.0, depth/20, depth/50]],  # 0.0 = auto (depth/10)
    :solver_tol => [(solver_tol=t,) for t ∈ [1e-4, 1e-6, 1e-8, 1e-10]],
    # combined coarse→fine configs: traces the speed-vs-accuracy frontier, so
    # we can ask "at BELLHOP-like run time, how much accuracy does RaySolver
    # lose?" (matched-time analysis)
    :frontier => [
      (nbeams=100, ds=depth/2, solver_tol=1e-3, min_amplitude=1e-3),
      (nbeams=250, ds=depth/2, solver_tol=1e-4, min_amplitude=1e-4),
      (nbeams=500, ds=depth/5, solver_tol=1e-4),
      (nbeams=1000, ds=depth/10, solver_tol=1e-6),
      (nbeams=2000, ds=depth/20, solver_tol=1e-8),
      (nbeams=4000, ds=depth/20, solver_tol=1e-8),
    ],
  )
end

"Bellhop option grids."
bellhop_sweeps(scn) = Dict(
  :nbeams => [(nbeams=n,) for n ∈ [250, 500, 1000, 2000, 4000, 0]],  # 0 = auto (0.05° spacing ⇒ 3201)
  :beam_type => [(beam_type=t,) for t ∈ [:geometric, :gaussian]],
)

# Converged ("best") configurations used for the headline comparison; refined
# manually after inspecting sweep results.
"BellhopJL (native port) option grids — mirrors the Fortran Bellhop sweeps."
bellhopjl_sweeps(scn) = Dict(
  :nbeams => [(nbeams=n,) for n ∈ [250, 500, 1000, 2000, 4000, 0]],
  :beam_type => [(beam_type=t,) for t ∈ [:geometric, :gaussian]],
)

best_opts(scn, ::Val{:raysolver}) = (nbeams=4000,)
best_opts(scn, ::Val{:bellhop}) = (nbeams=4000, beam_type=:gaussian)
best_opts(scn, ::Val{:bellhopjl}) = (nbeams=4000, beam_type=:gaussian)
best_opts(scn, m::Symbol) = best_opts(scn, Val(m))

