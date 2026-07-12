# Extracts benchmark results from the version-keyed cache into tidy CSVs under
# results/tables/<variant>/, where <variant> names the source state (e.g.
# "baseline" = benchmark branch, "pr26" = benchmark-pr26 branch). Run *on the
# corresponding branch* so cache keys match:
#
#   julia --project=. -t auto analyze.jl baseline
#
# report.qmd renders from these CSVs (plus figures), so it does not recompute.

include(joinpath(@__DIR__, "scenarios.jl"))
include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sweeps.jl"))
using CairoMakie

variant = isempty(ARGS) ? "baseline" : ARGS[1]
outdir = joinpath(RESULTS_DIR, "tables", variant)
mkpath(outdir)
figdir = joinpath(RESULTS_DIR, "figs", variant)
mkpath(figdir)

best_rs = (nbeams=4000,)
best_bh = (nbeams=4000, beam_type=:gaussian)
best_bjl = (nbeams=4000, beam_type=:gaussian)

# ---- best-vs-best and referee metrics ----
rows = NamedTuple[]
for scn ∈ SCENARIOS, mode ∈ (:incoherent, :coherent)
  wm = water_mask(scn)
  ta = tl_for(scn, :raysolver, best_rs; mode)
  tb = tl_for(scn, :bellhop, best_bh; mode)
  tj = try tl_for(scn, :bellhopjl, best_bjl; mode) catch e; @warn "bellhopjl tl failed" scn.name mode; nothing end
  push!(rows, (; scenario=scn.name, mode, pair="raysolver-vs-bellhop", tl_metrics(ta, tb; watermask=wm)...))
  if tj !== nothing
    push!(rows, (; scenario=scn.name, mode, pair="bellhopjl-vs-bellhop", tl_metrics(tj, tb; watermask=wm)...))
    push!(rows, (; scenario=scn.name, mode, pair="raysolver-vs-bellhopjl", tl_metrics(ta, tj; watermask=wm)...))
  end
  tr = try referee_tl(scn; mode) catch; nothing end
  if tr !== nothing
    push!(rows, (; scenario=scn.name, mode, pair="raysolver-vs-referee", tl_metrics(ta, tr; watermask=wm)...))
    push!(rows, (; scenario=scn.name, mode, pair="bellhop-vs-referee", tl_metrics(tb, tr; watermask=wm)...))
    tj !== nothing && push!(rows, (; scenario=scn.name, mode, pair="bellhopjl-vs-referee", tl_metrics(tj, tr; watermask=wm)...))
  end
end
CSV.write(joinpath(outdir, "best_vs_best.csv"), DataFrame(rows))

# ---- sweeps (self-convergence + frontier) ----
sw = DataFrame[]
for scn ∈ SCENARIOS
  for (msym, grids) ∈ [(:raysolver, raysolver_sweeps(scn)), (:bellhop, bellhop_sweeps(scn)),
                       (:bellhopjl, bellhopjl_sweeps(scn))]
    for (param, optslist) ∈ grids, mode ∈ (:incoherent, :coherent)
      df = sweep(scn, msym, optslist; mode)
      df.scenario .= scn.name; df.param .= string(param); df.mode .= string(mode)
      push!(sw, df)
    end
  end
end
CSV.write(joinpath(outdir, "sweeps.csv"), vcat(sw...))

# ---- timing ----
rows = NamedTuple[]
for scn ∈ SCENARIOS, msym ∈ (:raysolver, :bellhop, :bellhopjl)
  for (cfg, opts) ∈ [("default", (;)), ("best", msym === :raysolver ? best_rs : msym === :bellhop ? best_bh : best_bjl)]
    opts = msym === :bellhopjl && cfg == "best" ? best_bjl : opts
    b = try bench(scn, msym, opts) catch e
      @warn "bench failed" scn.name msym cfg
      continue
    end
    push!(rows, (; scenario=scn.name, model=msym, cfg, nthreads=Threads.nthreads(), b...))
  end
end
CSV.write(joinpath(outdir, "bench.csv"), DataFrame(rows))

# ---- arrivals comparison at point receivers ----
rows = NamedTuple[]
for scn ∈ SCENARIOS, (i, rx) ∈ enumerate(scn.rxs)
  a = try arrivals_for(scn, :raysolver, (;))[i] catch e
    @warn "raysolver arrivals failed" scn.name i
    missing
  end
  a === missing && continue
  b = arrivals_for(scn, :bellhop, (;))[i]
  m = match_arrivals(a, b)
  push!(rows, (; scenario=scn.name, rx=i, x=rx.pos.x, z=rx.pos.z,
               n_raysolver=length(a), n_bellhop=length(b), n_matched=nrow(m.matched),
               med_dt_us=nrow(m.matched) > 0 ? median(abs.(m.matched.dt_us)) : NaN,
               med_damp_db=nrow(m.matched) > 0 ? median(abs.(m.matched.damp_db)) : NaN))
end
CSV.write(joinpath(outdir, "arrivals.csv"), DataFrame(rows))

# ---- TL maps ----
for scn ∈ SCENARIOS, mode ∈ (:incoherent, :coherent)
  ta = tl_for(scn, :raysolver, best_rs; mode)
  tb = tl_for(scn, :bellhop, best_bh; mode)
  xs = collect(scn.grid.xrange) ./ 1000; zs = collect(scn.grid.zrange)
  lo = quantile(vec(filter(isfinite, ta)), 0.05)
  fig = Figure(size=(1000, 750))
  for (i, (t, lbl)) ∈ enumerate([(ta, "RaySolver"), (tb, "BELLHOP")])
    ax = Axis(fig[i,1], title="$lbl — $(scn.name) ($mode)", xlabel="range (km)", ylabel="depth (m)")
    hm = heatmap!(ax, xs, zs, t; colormap=Reverse(:jet), colorrange=(lo, lo+60))
    i == 1 && Colorbar(fig[1:2,2], hm, label="TL (dB)")
  end
  ax = Axis(fig[3,1], title="ΔTL (RaySolver − BELLHOP)", xlabel="range (km)", ylabel="depth (m)")
  hm = heatmap!(ax, xs, zs, ta .- tb; colormap=:balance, colorrange=(-10,10))
  Colorbar(fig[3,2], hm, label="ΔTL (dB)")
  save(joinpath(figdir, "$(scn.name)_$mode.png"), fig)
end

println("variant=$variant  version=$(version_key())  → $outdir")
