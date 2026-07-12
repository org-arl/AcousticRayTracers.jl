# Driver: runs all convergence sweeps, best-vs-best comparisons and timing
# benchmarks, populating benchmarks/results/cache. The report.qmd then renders
# from the cache. Run with:
#
#   julia --project=benchmarks -t auto benchmarks/run_benchmarks.jl [scenario ...]
#
# Passing scenario names limits the run. FORCE_RERUN=1 recomputes everything.

include(joinpath(@__DIR__, "scenarios.jl"))
include(joinpath(@__DIR__, "harness.jl"))

include(joinpath(@__DIR__, "sweeps.jl"))
# ---------------------------------------------------------------------------

function run_scenario(scn; modes=(:incoherent, :coherent))
  @info "=== scenario: $(scn.name)" scn.description
  for (msym, grids) ∈ [(:raysolver, raysolver_sweeps(scn)), (:bellhop, bellhop_sweeps(scn)),
                       (:bellhopjl, bellhopjl_sweeps(scn))]
    for (param, optslist) ∈ grids
      for mode ∈ modes
        @info "sweep" msym param mode
        df = sweep(scn, msym, optslist; mode)
        show(df, allrows=true, allcols=true); println()
      end
    end
  end
  # best-vs-best TL + arrivals + referee
  for mode ∈ modes
    ta = tl_for(scn, :raysolver, best_opts(scn, :raysolver); mode)
    tb = tl_for(scn, :bellhop, best_opts(scn, :bellhop); mode)
    @info "best-vs-best $(scn.name) $mode" tl_metrics(ta, tb)
    referee_tl(scn; mode)  # cache it if available
  end
  arrivals_for(scn, :raysolver, (;))
  arrivals_for(scn, :bellhop, (;))
  arrivals_for(scn, :bellhopjl, (;))
  # timing at defaults and at best
  for msym ∈ (:raysolver, :bellhop, :bellhopjl), opts ∈ [(;), best_opts(scn, msym)]
    @info "bench" msym opts bench(scn, msym, opts)
  end
end

names = isempty(ARGS) ? [s.name for s ∈ SCENARIOS] : ARGS
for n ∈ names
  try
    run_scenario(scenario(n))
  catch err
    @error "scenario $n failed" err
  end
end
@info "done" version_key()
