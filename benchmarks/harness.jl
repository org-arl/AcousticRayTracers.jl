# Shared harness for the RaySolver vs BELLHOP benchmark.
#
# Provides: model construction from an options NamedTuple, TL computation,
# comparison metrics, arrival matching, convergence-sweep and timing drivers,
# and disk caching keyed by the source versions in use (so results
# auto-invalidate when AcousticRayTracers or UnderwaterAcoustics change).

using UnderwaterAcoustics
using AcousticRayTracers
using AcousticsToolbox
using BellhopJL
using Statistics, Printf, Serialization
using DataFrames, CSV
using BenchmarkTools

const RESULTS_DIR = joinpath(@__DIR__, "results")
const CACHE_DIR = joinpath(RESULTS_DIR, "cache")

# ---------------------------------------------------------------------------
# Versioning / caching

"Short git SHA of the AcousticRayTracers source being benchmarked."
art_sha() = readchomp(`git -C $(dirname(@__DIR__)) rev-parse --short HEAD`)

"Identifier for the UnderwaterAcoustics version in the manifest."
function uwa_rev()
  m = read(joinpath(@__DIR__, "Manifest.toml"), String)
  b = match(r"\[\[deps\.UnderwaterAcoustics\]\][^[]*"s, m).match
  h = match(r"git-tree-sha1 = \"([0-9a-f]+)\"", b)
  h === nothing ? "reg" : h[1][1:8]
end

version_key() = "$(art_sha())-$(uwa_rev())"

"Short git SHA of the BellhopJL source (native Julia Bellhop port)."
bjl_sha() = try
  readchomp(`git -C /Users/mandar/Projects/BellhopJL.jl rev-parse --short HEAD`)
catch
  "none"
end

"""
    cached(f, parts...; force=false)

Run `f()` and cache its result under a key derived from `parts` and the
current source versions. Set `force=true` (or ENV["FORCE_RERUN"]="1") to
recompute.
"""
function cached(f, parts...; force=false)
  mkpath(CACHE_DIR)
  # BellhopJL results are keyed additionally by its own source SHA
  vkey = any(p -> p === :bellhopjl, parts) ? "$(version_key())-$(bjl_sha())" : version_key()
  key = string(hash((vkey, parts...)); base=16)
  path = joinpath(CACHE_DIR, "$key.jls")
  if !force && get(ENV, "FORCE_RERUN", "0") != "1" && isfile(path)
    return deserialize(path)
  end
  x = f()
  serialize(path, x)
  x
end

# ---------------------------------------------------------------------------
# Model construction

"""
    make_model(:raysolver | :bellhop, env; opts...)

Construct a propagation model from an options NamedTuple. Unknown options for
the given model raise an error (so sweeps can't silently no-op).
"""
make_model(m::Symbol, env; opts...) = make_model(Val(m), env; opts...)
make_model(::Val{:raysolver}, env; opts...) = RaySolver(env; opts...)
make_model(::Val{:bellhop}, env; opts...) = Bellhop(env; opts...)
make_model(::Val{:bellhopjl}, env; opts...) = BellhopModel(env; opts...)

# ---------------------------------------------------------------------------
# TL computation and metrics

"Transmission loss (dB) over the scenario grid."
function compute_tl(model, tx, grid; mode=:incoherent)
  x = acoustic_field(model, tx, grid; mode)
  -20.0 .* log10.(abs.(x) .+ eps())
end

"Cached TL for (scenario, model symbol, options, mode)."
function tl_for(scn, msym::Symbol, opts::NamedTuple; mode=:incoherent, force=false)
  cached(scn.name, msym, opts, mode; force) do
    model = make_model(msym, scn.env; opts...)
    compute_tl(model, scn.tx, scn.grid; mode)
  end
end

"""
    tl_metrics(tl_a, tl_b; maxtl=100.0)

Pointwise |ΔTL| stats over grid points where both models predict TL below
`maxtl` dB (masks out shadow zones / numerical floors).
"""
function tl_metrics(tl_a::AbstractMatrix, tl_b::AbstractMatrix; maxtl=100.0, watermask=trues(size(tl_a)))
  mask = isfinite.(tl_a) .& isfinite.(tl_b) .& (tl_a .< maxtl) .& (tl_b .< maxtl) .& watermask
  d = abs.(tl_a[mask] .- tl_b[mask])
  isempty(d) && return (median=NaN, p90=NaN, max=NaN, within3dB=NaN, npts=0)
  (median=median(d), p90=quantile(d, 0.9), max=maximum(d),
   within3dB=mean(d .< 3.0), npts=length(d))
end

"""
    water_mask(scn)

BitMatrix over the scenario grid marking receivers inside the water column
(above the bathymetry at their range). BELLHOP deposits finite field into
below-seafloor receivers (halfspace), which would otherwise pollute metrics in
range-dependent-bathymetry scenarios.
"""
function water_mask(scn)
  [ -value(scn.env.bathymetry, (x, 0.0, 0.0)) < z
    for x ∈ scn.grid.xrange, z ∈ scn.grid.zrange ]
end

# ---------------------------------------------------------------------------
# Arrivals comparison

"Cached arrivals at each scenario point receiver."
function arrivals_for(scn, msym::Symbol, opts::NamedTuple; force=false)
  cached(scn.name, msym, opts, :arrivals; force) do
    model = make_model(msym, scn.env; opts...)
    [arrivals(model, scn.tx, rx; paths=false) for rx ∈ scn.rxs]
  end
end

"""
    match_arrivals(a, b; ttol=1e-3)

Match arrival lists by (surface bounces, bottom bounces, nearest travel time).
Returns a DataFrame of matched pairs (Δt in µs, Δamp in dB, Δphase in deg) and
counts of unmatched arrivals on each side.
"""
function match_arrivals(a, b; ttol=1e-3)
  used = falses(length(b))
  rows = NamedTuple[]
  for x ∈ a
    cand = [j for j ∈ eachindex(b) if !used[j] && b[j].ns == x.ns && b[j].nb == x.nb &&
            abs(b[j].t - x.t) < ttol]
    isempty(cand) && continue
    j = argmin(j -> abs(b[j].t - x.t), cand)
    used[j] = true
    y = b[j]
    push!(rows, (ns=x.ns, nb=x.nb, t=x.t,
                 dt_us=(y.t - x.t) * 1e6,
                 damp_db=20log10(abs(y.ϕ) + eps()) - 20log10(abs(x.ϕ) + eps()),
                 dphase_deg=rad2deg(angle(y.ϕ / x.ϕ))))
  end
  (matched=DataFrame(rows), only_a=length(a) - length(rows), only_b=count(!, used))
end

"""
    ray_paths_for(scn, msym, rx; opts=(;), force=false)

Cached eigenray arrivals *with paths* at a given receiver — used to attribute
significant TL differences: if matched eigenrays have the same geometry and
travel time but different amplitudes, the discrepancy is in amplitude/beam
computation; if the paths themselves diverge, it is in the ray trace.
"""
function ray_paths_for(scn, msym::Symbol, rx; opts=(;), force=false)
  cached(scn.name, msym, opts, :raypaths, (rx.pos.x, rx.pos.z); force) do
    model = make_model(msym, scn.env; opts...)
    arrivals(model, scn.tx, rx; paths=true)
  end
end

"""
    ray_diagnostics(scn, rx; opts_a=(;), opts_b=(;))

Match eigenrays between the two models at a receiver and report, per matched
pair: launch/arrival angles, bounce counts, Δ travel time, Δ amplitude, and the
max spatial divergence between the two paths (interpolated onto common ranges).
Separates ray-trace errors (geometry/time) from field errors (amplitude).
"""
function ray_diagnostics(scn, rx; opts_a=(;), opts_b=(;))
  a = ray_paths_for(scn, :raysolver, rx; opts=opts_a)
  b = ray_paths_for(scn, :bellhop, rx; opts=opts_b)
  used = falses(length(b))
  rows = NamedTuple[]
  for x ∈ a
    cand = [j for j ∈ eachindex(b) if !used[j] && b[j].ns == x.ns && b[j].nb == x.nb]
    isempty(cand) && continue
    j = argmin(j -> abs(b[j].t - x.t), cand)
    used[j] = true
    y = b[j]
    # max |Δz| between paths sampled at common ranges
    dz = NaN
    if x.path !== missing && y.path !== missing && length(x.path) > 1 && length(y.path) > 1
      xr = [p.x for p ∈ x.path]; xz = [p.z for p ∈ x.path]
      yr = [p.x for p ∈ y.path]; yz = [p.z for p ∈ y.path]
      rs = range(max(minimum(xr), minimum(yr)), min(maximum(xr), maximum(yr)); length=200)
      interp(r, rr, zz) = begin
        i = clamp(searchsortedlast(rr, r), 1, length(rr)-1)
        α = clamp((r - rr[i]) / (rr[i+1] - rr[i] + eps()), 0.0, 1.0)
        zz[i] + α * (zz[i+1] - zz[i])
      end
      issorted(xr) && issorted(yr) &&
        (dz = maximum(abs(interp(r, xr, xz) - interp(r, yr, yz)) for r ∈ rs))
    end
    push!(rows, (ns=x.ns, nb=x.nb,
                 θs_a=rad2deg(x.θₛ), θs_b=rad2deg(y.θₛ),
                 t_a_ms=x.t*1e3, dt_us=(y.t - x.t)*1e6,
                 amp_a_db=20log10(abs(x.ϕ)+eps()), damp_db=20log10(abs(y.ϕ)+eps()) - 20log10(abs(x.ϕ)+eps()),
                 max_dz_m=dz))
  end
  (matched=DataFrame(rows),
   only_raysolver=[(ns=x.ns, nb=x.nb, t=x.t, θₛ=rad2deg(x.θₛ)) for x ∈ a if !any(r -> r.ns == x.ns && r.nb == x.nb && abs(r.t_a_ms - x.t*1e3) < 1e-9, rows)],
   only_bellhop=[(ns=b[j].ns, nb=b[j].nb, t=b[j].t, θₛ=rad2deg(b[j].θₛ)) for j ∈ eachindex(b) if !used[j]])
end

"Overlay eigenray paths from both models at a receiver (RaySolver solid blue, BELLHOP dashed red)."
function ray_overlay_figure(scn, rx; opts_a=(;), opts_b=(;))
  a = ray_paths_for(scn, :raysolver, rx; opts=opts_a)
  b = ray_paths_for(scn, :bellhop, rx; opts=opts_b)
  fig = Figure(size=(1000, 400))
  ax = Axis(fig[1,1], title="eigenrays — $(scn.name) at ($(rx.pos.x), $(rx.pos.z))",
            xlabel="range (m)", ylabel="depth (m)")
  for x ∈ a
    x.path === missing && continue
    lines!(ax, [p.x for p ∈ x.path], [p.z for p ∈ x.path]; color=(:blue, 0.5))
  end
  for y ∈ b
    y.path === missing && continue
    lines!(ax, [p.x for p ∈ y.path], [p.z for p ∈ y.path]; color=(:red, 0.6), linestyle=:dash)
  end
  fig
end

# ---------------------------------------------------------------------------
# Convergence sweep

"""
    sweep(scn, msym, optslist; mode=:incoherent)

Run the model at each options NamedTuple in `optslist`; compare every result
against the last entry (assumed the finest/most converged). Also records a
single-run wall time (coarse; use `bench` for careful timing). Returns a
DataFrame with one row per configuration.
"""
function sweep(scn, msym::Symbol, optslist; mode=:incoherent, maxtl=100.0)
  ref = tl_for(scn, msym, optslist[end]; mode)
  rows = NamedTuple[]
  for opts ∈ optslist
    t = cached(scn.name, msym, opts, mode, :walltime) do
      model = make_model(msym, scn.env; opts...)
      compute_tl(model, scn.tx, scn.grid; mode)  # warm-up
      @elapsed compute_tl(model, scn.tx, scn.grid; mode)
    end
    m = tl_metrics(tl_for(scn, msym, opts; mode), ref; maxtl)
    push!(rows, (; model=msym, opts=repr(opts), m..., walltime_s=t))
  end
  DataFrame(rows)
end

# ---------------------------------------------------------------------------
# Careful timing (for the final, converged configurations)

"""
    bench(scn, msym, opts; mode=:incoherent, seconds=20)

BenchmarkTools timing of the full acoustic_field grid computation, and of
single-receiver arrivals. Returns (field_ms_min, field_ms_median, arr_ms_min,
arr_ms_median).
"""
function bench(scn, msym::Symbol, opts::NamedTuple; mode=:incoherent, seconds=20)
  cached(scn.name, msym, opts, mode, :bench, Threads.nthreads()) do
    model = make_model(msym, scn.env; opts...)
    rx1 = scn.rxs[1]
    bf = @benchmark acoustic_field($model, $(scn.tx), $(scn.grid); mode=$mode) seconds=seconds
    ba = @benchmark arrivals($model, $(scn.tx), $rx1; paths=false) seconds=seconds
    (field_ms_min=minimum(bf).time / 1e6, field_ms_median=median(bf).time / 1e6,
     arr_ms_min=minimum(ba).time / 1e6, arr_ms_median=median(ba).time / 1e6)
  end
end

# ---------------------------------------------------------------------------
# Referees

"""
TL from the referee model for a scenario (or nothing). Note: full-field models
(modes/wavenumber-integration) are inherently coherent — referee comparisons
should use `mode=:coherent`, or compare range-averaged TL for an
incoherent-like smoothing.
"""
function referee_tl(scn; mode=:coherent)
  scn.referee === nothing && return nothing
  cached(scn.name, scn.referee, (;), mode, :referee) do
    model = scn.referee === :pekeris ? PekerisModeSolver(scn.env) :
            scn.referee === :kraken ? Kraken(scn.env) :
            scn.referee === :lloyd ? PekerisRayTracer(scn.env) :
            error("unknown referee $(scn.referee)")
    compute_tl(model, scn.tx, scn.grid; mode)
  end
end
