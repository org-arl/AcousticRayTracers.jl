import LinearAlgebra: norm, dot
import OrdinaryDiffEq: ODEProblem, VectorContinuousCallback, CallbackSet, Tsit5, solve
import OrdinaryDiffEqRosenbrock: Rosenbrock23
import DiffEqCallbacks: StepsizeLimiter
import NonlinearSolve: IntervalNonlinearProblem, NonlinearProblem
import SciMLBase: successful_retcode, terminate!, remake
import ForwardDiff
import StaticArrays: SA, SVector

# scatterers from the environment are wrapped into a ScattererSet at solver
# construction; geometry queries go through the UnderwaterAcoustics scatterer API
# (distance, boundary_projection, reflection_coef)
struct ScattererSet{V}
  items::V
end

Base.@kwdef struct RaySolver{T1,T2,T3} <: AbstractRayPropagationModel
  env::T1
  nbeams::Int = 0
  min_angle::Float64 = -deg2rad(80)
  max_angle::Float64 = +deg2rad(80)
  ds::Float64 = 0.0
  atol::Float64 = 1e-4
  rugosity::Float64 = 1.5
  min_amplitude::Float64 = 1e-6
  solver::T2 = nothing
  solver_tol::Float64 = 1e-8
  backscatter::Bool = false
  rmin::Float64 = 0.0
  rmax::Float64 = 0.0
  scatterers::T3 = nothing
  ntasks::Int = 0
  function RaySolver(env, nbeams, min_angle, max_angle, ds, atol, rugosity, min_amplitude, solver, solver_tol, backscatter, rmin, rmax, scatterers, ntasks)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ min_angle ≤ π/2 || error("min_angle should be between -π/2 and π/2")
    -π/2 ≤ max_angle ≤ π/2 || error("max_angle should be between -π/2 and π/2")
    min_angle < max_angle || error("max_angle should be more than min_angle")
    rmin ≤ 0 || error("rmin should be non-positive")
    rmax ≥ 0 || error("rmax should be non-negative")
    solver = something(solver, is_isovelocity(env) ? Tsit5() : Rosenbrock23())
    scatterers === nothing && (scatterers = _scatterer_set(env))
    new{typeof(env),typeof(solver),typeof(scatterers)}(env, nbeams, min_angle, max_angle, ds, atol, rugosity, min_amplitude, solver, solver_tol, backscatter, rmin, rmax, scatterers, ntasks)
  end
end

"""
    RaySolver(env; kwargs...)

Julia implementation of a ray/Gaussian beam propagation model. The model supports
complex environments, but retains differentiability.

Supported keyword arguments:
- `nbeams`: number of beams to use (default: 0, auto)
- `min_angle`: minimum beam angle (default: -80°)
- `max_angle`: maximum beam angle (default: 80°)
- `ds`: nominal spacing between ray points (default: 1/10 water depth)
- `atol`: absolute position tolerance (default: 0.0001 m)
- `rugosity`: TODO (default: 1.5)
- `min_amplitude`: minimum ray amplitude to track (default: 1e-6); rays weaker
  than this (including absorption) stop being traced, which also culls weak
  multipath from `arrivals` — set to `0.0` to keep all arrivals (e.g. for
  impulse-response work)
- `solver`: differential equation solver (default: nothing, auto)
- `solver_tol`: solver tolerance (default: 1e-8)
- `backscatter`: continue tracing rays that turn back toward the source (default: false)
- `rmin`: left edge of modeling domain in m, ≤ 0 (default: 0; only meaningful with backscatter)
- `rmax`: right edge of modeling domain in m (default: 0, auto = query range)
- `ntasks`: number of parallel tasks to use (default: 0, auto = number of threads;
  1 = serial, useful for benchmarking and for AD tools that dislike threads)

If the environment contains scatterers (`env.scatterers`), rays reflect off them
using the scatterer's boundary condition. With `backscatter` enabled, rays are
traced until they exit the modeling domain `r ∈ [rmin, rmax]` or become too weak
(`min_amplitude`, including absorption); backscattered arrivals are reported with
|arrival angle| > π/2. Features outside the modeling domain are invisible; set
`rmax` beyond the query range to capture echoes from behind the receivers.
If `rmin < 0`, a second fan of rays with the same elevation angles is launched
in the -x direction, so that features to the left of the transmitter can also
produce echoes (such arrivals have |launch angle| > π/2).
"""
RaySolver(env; kwargs...) = RaySolver(; env=env, kwargs...)

Base.show(io::IO, pm::RaySolver) = print(io, "RaySolver(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::RaySolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true, ztol=0.1)
  _check2d([tx], [rx])
  pm.backscatter && return _arrivals_backscatter(pm, tx, rx; paths, ztol)
  p2 = location(rx)
  p1 = location(tx)
  R = abs(p2.x - p1.x)
  nbeams = pm.nbeams
  if nbeams == 0
    h = maximum(pm.env.bathymetry)
    nbeams = ceil(Int, 16 * (pm.max_angle - pm.min_angle) / atan(h, R))
  end
  ds = pm.ds ≤ 0 ? minimum(pm.env.bathymetry) / 10 : pm.ds
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  # eigenrays are found from receiver-range crossings grouped by bounce
  # signature (_fan_search): with range-dependent bathymetry, adjacent fan rays
  # can follow different bounce histories (e.g. hit or clear an obstacle), so
  # the raw depth error Δz(θ) is discontinuous and brackets across a signature
  # change are meaningless; within one signature Δz(θ) is continuous. Using
  # crossings (rather than ray endpoints) also finds arrivals on rays that
  # refract back down-range after passing the receiver.
  rdom = p2.x + 0.1   # trace slightly past rx so the range crossing is recorded before termination
  T1 = promote_type(env_type(pm.env), eltype(location(tx)), typeof(p2.x), typeof(p2.z))
  T2 = typeof(RayArrival(zero(T1), zero(Complex{T1}), 0, 0, zero(T1), zero(T1),
    Array{@NamedTuple{x::T1,y::T1,z::T1}}(undef, 0)))
  erays = _fan_search(T2, pm, tx, collect(T1, θ), T1(rdom), T1(p2.x), T1(p2.z), ds, ztol, paths; edges=true)
  sort!(erays; by=Base.Fix2(getfield, :t))
end

function UnderwaterAcoustics.acoustic_field(pm::RaySolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; mode=:coherent)
  arr = arrivals(pm, tx, rx)
  length(arr) == 0 && return zero(UnderwaterAcoustics._phasortype(eltype(arr)))
  if mode === :incoherent
    complex(√sum(a -> abs2(a.ϕ), arr)) * db2amp(spl(tx))
  elseif mode === :coherent
    f = frequency(tx)
    sum(a -> a.ϕ * cispi(2f * a.t), arr) * db2amp(spl(tx))
  else
    error("Unknown mode :$mode")
  end
end

# Gaussian beam implementation primarily based on ideas from:
# [COA – Computational Ocean Acoustics, 2nd ed., ch. 3)]
function UnderwaterAcoustics.acoustic_field(pm::RaySolver, tx::AbstractAcousticSource, rxs::AcousticReceiverGrid2D; mode=:coherent)
  _check2d([tx], rxs)
  mode ∈ (:coherent, :incoherent) || error("Unknown mode :" * string(mode))
  if pm.nbeams > 0
    nbeams = pm.nbeams
  else
    p1 = location(tx)
    R = abs(maximum(rxs.xrange) - p1[1])
    Δz = abs(Float64(rxs.zrange.step))
    # cap was 1000, but benchmark convergence sweeps vs BELLHOP show deep-water
    # grids need ~4000 beams (~0.04°/beam over a ±80° fan) to self-converge;
    # BELLHOP's own auto uses 0.05° spacing
    nbeams = clamp(ceil(Int, 8 * (pm.max_angle - pm.min_angle) * R / Δz), 100, 4000)
  end
  ds = pm.ds ≤ 0 ? minimum(pm.env.bathymetry) / 10 : pm.ds
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  δθ = Float64(θ.step)
  f = frequency(tx)
  h = maximum(pm.env.bathymetry)
  mid_temp = (minimum(pm.env.temperature) + maximum(pm.env.temperature)) / 2
  c₀ = value(pm.env.soundspeed, location(tx))
  rmax = maximum(rxs.xrange) + 0.1
  pm.backscatter && pm.rmax > 0 && (rmax = pm.rmax)
  γ = absorption(f, 1, pm.env.salinity, mid_temp, h/2)  # nominal absorption per meter
  log_γ = log(γ)
  T = promote_type(env_type(pm.env), eltype(location(tx)), typeof(f), ComplexF64)
  θs = pm.backscatter && pm.rmin < 0 ? vcat(collect(θ), rem2pi.(π .- θ, RoundNearest)) : collect(θ)
  # deposit the Gaussian beam contributions of one ray onto a field array;
  # each task accumulates into its own private array (no locks), and the
  # per-task arrays are summed at the end
  function beam!(afld, θ₀)
    ray, aux_info = _trace(pm, tx, θ₀, rmax, ds; aux=true)
    for i ∈ 1:length(ray.path)-1
      pos1 = SA[ray.path[i].x, ray.path[i].z]           # start of ray segment
      vec12 = SA[ray.path[i+1].x, ray.path[i+1].z] - pos1   # vector to end of ray segment
      vec12_mag2 = dot(vec12, vec12)
      s1, q1, t1, A1 = aux_info[i]
      s2, q2, t2, _ = aux_info[i+1]
      γₛ = exp(log_γ * (s1 + s2) / 2)                   # nominal absorption integrated
      cₛ = let p = pos1 + vec12 / 2                     # midpoint of ray segment
        value(pm.env.soundspeed, (p[1], 0.0, p[2]))     # nominal soundspeed
      end
      C1 = A1 * γₛ / √π                                 # √π scaling for Gaussian beams
      mode === :incoherent && (C1 *= π/2)               # scaling for incoherent summation
      C2 = δθ * abs(cos(θ₀)) * cₛ / c₀
      xi, zi = _findall_rxgrid2d(rxs,                   # subset of rx in neighborhood
        ray.path[i].x, ray.path[i+1].x,                 #  of ray segment
        ray.path[i].z, ray.path[i+1].z,
        # neighborhood must cover the clamped minimum beam width (≤ πλ) too
        max(4 * max(abs(q1), abs(q2)) * δθ, 4π * cₛ / f)
      )
      for j ∈ xi                                        # loop over the subset
        zbot = -value(pm.env.bathymetry, (rxs[j,1].pos.x, 0.0, 0.0))
        ztop = value(pm.env.altimetry, (rxs[j,1].pos.x, 0.0, 0.0))
        for k ∈ zi
          rx = rxs[j,k].pos
          zbot ≤ rx.z ≤ ztop || continue                # no deposit outside the water column
          rxpos = SA[rx.x, rx.z]
          α = dot(rxpos - pos1, vec12) / vec12_mag2
          0 ≤ α < 1 || continue                         # rx outside of ray segment
          q = q1 + α * (q2 - q1)                        # spreading factor at cpa
          t = t1 + α * (t2 - t1)                        # time at cpa
          W = abs(q * δθ)                               # beam width at cpa [COA (3.74)]
          # in incoherent mode, clamp beam width from below near caustics (q → 0
          # makes the beam infinitesimally narrow and its amplitude blow up);
          # same rule (and units quirk) as BELLHOP's InfluenceGeoGaussianCart:
          # min(0.2 f t, π λ). Not applied in coherent mode, where beams wider
          # than the multipath image separation would blur interference fringes.
          mode === :incoherent && (W = max(W, min(0.2 * f * t, π * cₛ / f)))
          cpa = pos1 + α * vec12                        # closest point of approach
          cpa[1] > 0 || continue                        # no deposit at/behind the r=0 axis
          n = norm(cpa - rxpos)                         # normal distance from ray to rx
          n < 4W || continue                            # rx too far from ray segment
          A = C1 * sqrt(C2 / (cpa[1] * W))              # [COA (3.76)]
          if mode === :coherent
            P = A * exp(-(n / W)^2) * cispi(2f * t)     # [COA (3.72), COA (3.75)]
          else
            P = complex(abs2(A * exp(-(n / W)^2)), 0.0)
          end
          afld[j,k] += P
        end
      end
    end
    afld
  end
  nt = min(_ntasks(pm), length(θs))
  afld = zeros(T, size(rxs,1), size(rxs,2))
  if nt ≤ 1
    foreach(θ₀ -> beam!(afld, θ₀), θs)
  else
    # beams are interleaved across tasks for load balancing (work per beam is
    # uneven); each task returns its private field array
    tasks = map(1:nt) do c
      Threads.@spawn begin
        buf = zeros(T, size(rxs,1), size(rxs,2))
        for θ₀ ∈ @view θs[c:nt:end]
          beam!(buf, θ₀)
        end
        buf
      end
    end
    for t ∈ tasks
      afld .+= fetch(t)::Matrix{T}
    end
  end
  mode === :incoherent && (afld .= sqrt.(afld))
  afld * db2amp(spl(tx))
end

### helper functions

_ordered(a, b) = a < b ? (a, b) : (b, a)

# effective number of parallel tasks (0 = auto)
_ntasks(pm::RaySolver) = pm.ntasks ≤ 0 ? Threads.nthreads() : pm.ntasks

# run f(i) for i ∈ 1:n over up to ntasks tasks; iterations are interleaved
# across tasks for load balancing, since work per iteration is typically very
# uneven; ntasks ≤ 1 gives a fully serial code path (AD-friendly)
function _tforeach(f, n::Int, ntasks::Int)
  nt = min(ntasks, n)
  if nt ≤ 1
    for i ∈ 1:n
      f(i)
    end
  else
    @sync for c ∈ 1:nt
      Threads.@spawn for i ∈ c:nt:n
        f(i)
      end
    end
  end
  nothing
end

# serial-aware version of tmap
_tmap(f, x, ntasks) = ntasks ≤ 1 ? map(f, x) : tmap(f, x)

function _isnearzero(a, b, c)
  (isnan(a) || isnan(b) || isnan(c)) && return false
  sign(a) == sign(b) == sign(c) || return false
  abs(a) < abs(b) && return false
  abs(c) < abs(b) && return false
  return true
end

function _check2d(tx, rx)
  all(location(tx1).x == 0.0 for tx1 ∈ tx) || error("RaySolver requires transmitters at (0, 0, z)")
  all(location(tx1).y == 0.0 for tx1 ∈ tx) || error("RaySolver requires transmitters in the x-z plane")
  all(location(rx1).x >= 0.0 for rx1 ∈ rx) || error("RaySolver requires receivers to be in the +x halfspace")
  all(location(rx1).y == 0.0 for rx1 ∈ rx) || error("RaySolver requires receivers in the x-z plane")
end

function _findall_rxgrid2d(rxs, r1, r2, z1, z2, tol)
  rmin, rmax = _ordered(r1, r2)
  zmin, zmax = _ordered(z1, z2)
  rmin -= tol
  rmax += tol
  zmin -= tol
  zmax += tol
  x0 = first(rxs.xrange)
  dx = Float64(rxs.xrange.step)
  i = max(ceil(Int, (rmin - x0) / dx) + 1, 1)
  j = min(floor(Int, (rmax - x0) / dx) + 1, length(rxs.xrange))
  z0 = first(rxs.zrange)
  dz = Float64(rxs.zrange.step)
  k = max(ceil(Int, (zmin - z0) / dz) + 1, 1)
  l = min(floor(Int, (zmax - z0) / dz) + 1, length(rxs.zrange))
  max(min(i,j),1):min(max(i,j),length(rxs.xrange)), max(min(k,l),1):min(max(k,l),length(rxs.zrange))
end

### scatterer support

# wrap env.scatterers (UnderwaterAcoustics scatterer API) into a ScattererSet;
# nothing if the environment has no scatterers
function _scatterer_set(env)
  s = hasproperty(env, :scatterers) ? env.scatterers : ()
  length(s) == 0 ? nothing : ScattererSet(s)
end

# element type carried by scatterer geometry/boundary (may be a dual number)
_scat_eltype(::Nothing) = Bool
_scat_eltype(ss::ScattererSet) = mapreduce(_obj_eltype, promote_type, ss.items; init=Bool)
_obj_eltype(x::Number) = typeof(real(x))
_obj_eltype(x::Union{AbstractString,Function,Symbol}) = Bool
_obj_eltype(x) = nfields(x) == 0 ? Bool : mapreduce(i -> _obj_eltype(getfield(x, i)), promote_type, 1:nfields(x); init=Bool)

# minimum signed distance from (r, z) to any scatterer boundary (1-Lipschitz)
_mindist(::Nothing, r, z) = nothing
function _mindist(ss::ScattererSet, r, z)
  pos = (x=r, y=zero(r), z=z)
  minimum(s -> UnderwaterAcoustics.distance(s.shape, pos), ss.items)
end

# index of the scatterer nearest to (r, z)
function _nearest(ss::ScattererSet, r, z)
  pos = (x=r, y=zero(r), z=z)
  jbest = 1
  dbest = UnderwaterAcoustics.distance(ss.items[1].shape, pos)
  for j ∈ 2:length(ss.items)
    d = UnderwaterAcoustics.distance(ss.items[j].shape, pos)
    d < dbest && (dbest = d; jbest = j)
  end
  jbest
end

### ray tracing core

function _∂u((r, z, ξ, ζ, t, p, q), (c, ∂c, ∂²c), s)
  # implementation based on [COA (3.161-164, 3.58-63)]
  # assumes range-independent soundspeed
  cᵥ = c(z)
  cᵥ² = cᵥ * cᵥ
  c̄ = ∂²c(z) * ξ * ξ
  SA[cᵥ * ξ, cᵥ * ζ, 0, -∂c(z) / cᵥ², 1 / cᵥ, -c̄ * q, cᵥ * p]
end

# events: 1 = surface, 2 = bottom, 3 = max range, 4 = ray turned back (legacy) or
# exited past the source (backscatter), 5 = scatterer hit, 6 = receiver range
# crossing (recorded, not terminal), 6+k = crossing of the k-th soundspeed-
# gradient discontinuity (transmission correction applied, not terminal),
# 6+nk+m = crossing of the m-th boundary range-knot (no-op, not terminal: it
# forces an integration step exactly at the knot so that a boundary vertex —
# e.g. a seamount peak — cannot be straddled by a single step whose endpoints
# are both in the water, which would let the ray tunnel through the feature).
# Inert slots hold one(u[1]) so they never fire and preserve the element type
# (dual numbers).
function _check_ray!(out, u, s, integrator, a, b, rmax, rmin, backscatter, ss, δ, guard, r_rx, knots, rknots)
  out[1] = value(a, (u[1], 0.0, 0.0)) - u[2]    # surface reflection
  out[2] = u[2] + value(b, (u[1], 0.0, 0.0))    # bottom reflection
  out[3] = rmax - u[1]                          # maximum range
  out[4] = backscatter ? u[1] - rmin + δ : u[3] # exit domain left / ray turned back
  out[5] = ss === nothing || s < guard ? one(u[1]) : _mindist(ss, u[1], u[2]) - δ
  out[6] = isnan(r_rx) ? one(u[1]) : u[1] - r_rx
  nk = 0
  if knots !== nothing
    nk = length(knots.list)
    for k ∈ 1:nk
      out[6+k] = u[2] - knots.list[k][1]
    end
  end
  for m ∈ eachindex(rknots)
    out[6+nk+m] = u[1] - rknots[m]
  end
end

# handle a fired ray event: event 6 (receiver-range crossing) is recorded and
# integration continues; events 7+ (soundspeed-gradient discontinuity crossings)
# apply the transmission correction to p and continue; any other event terminates
# the segment and its index is reported via evref. Depending on the DiffEqBase
# version, ndx is either the event index or a vector of simultaneous event flags
# (0 = not fired, ±1 = fired). If a terminal event fires simultaneously with a
# knot crossing (e.g. a bounce off a boundary at knot depth), the ray reflects
# rather than transmits, so no correction is applied.
function _ray_event!(i, ndx::Integer, crossbuf, evref, knots, a, b)
  nk = knots === nothing ? 0 : length(knots.list)
  if ndx == 6
    push!(crossbuf, (i.t, i.u))
  elseif ndx > 6 + nk
    # boundary range-knot crossing: no action needed, the step landing here is
    # the point (see _check_ray!)
  elseif ndx > 6
    _knot_jump!(i, knots.list[ndx-6], a, b, evref)
  else
    evref[] = Int(ndx)
    terminate!(i)
  end
end

function _ray_event!(i, ndx::AbstractVector, crossbuf, evref, knots, a, b)
  nk = knots === nothing ? 0 : length(knots.list)
  length(ndx) ≥ 6 && ndx[6] != 0 && push!(crossbuf, (i.t, i.u))
  for k ∈ 1:min(length(ndx), 5)
    if ndx[k] != 0
      evref[] = k
      terminate!(i)
      return
    end
  end
  for k ∈ 7:min(length(ndx), 6 + nk)
    ndx[k] != 0 && _knot_jump!(i, knots.list[k-6], a, b, evref)
  end
end

# depths at which the soundspeed gradient is discontinuous, for piecewise-linear
# (SampledField ... interp=Linear()) depth-only SSPs. The dynamic ray tracing
# equations integrate p using ∂²c/∂z², which contains a delta at such depths that
# the ODE integration cannot see; a transmission correction to p is applied there
# instead (see _knot_jump!). Returns nothing (no events needed) for smooth
# profiles. list entries are (z_k, g₋, g₊, c_k) with g₋/g₊ the SSP slopes
# below/above z_k (zero outside the sampled domain, per flat extrapolation) and
# c_k = c(z_k); κ bounds the ray curvature (max |∂c|/c) and cmax the soundspeed,
# both used by the knot step-size limiter. Knots at or outside the boundaries
# are excluded: rays cannot cross them.
function _gradient_discontinuities(env)
  ss = env.soundspeed
  ss isa SampledFieldZ || return nothing
  ss.interp isa UnderwaterAcoustics.Linear || return nothing
  zs = sort([Float64(z) for z ∈ ss.zrange])
  cs = [value(ss, z) for z ∈ zs]
  n = length(zs)
  slopes = [(cs[k+1] - cs[k]) / (zs[k+1] - zs[k]) for k ∈ 1:n-1]
  zmax = minimum(env.altimetry) - 1e-3
  zmin = -maximum(env.bathymetry) + 1e-3
  list = [(zs[k], k == 1 ? zero(slopes[1]) : slopes[k-1],
           k == n ? zero(slopes[1]) : slopes[k], cs[k]) for k ∈ 1:n]
  filter!(list) do (z, g₋, g₊, _)
    zmin < z < zmax && !(g₋ ≈ g₊)
  end
  isempty(list) && return nothing
  (; list, κ=maximum(abs, slopes) / minimum(cs), cmax=maximum(cs))
end

# ranges at which the water-column boundaries have slope discontinuities
# (interior sample points of piecewise-linear SampledField bathymetry or
# altimetry). Used as no-op ray events so that an integration step cannot
# straddle a boundary vertex (see _check_ray!). Returns an empty vector when
# both boundaries are smooth in range.
function _boundary_range_knots(env)
  rknots = Float64[]
  for fld ∈ (env.bathymetry, env.altimetry)
    fld isa SampledFieldX || continue
    fld.interp isa UnderwaterAcoustics.Linear || continue
    xs = sort([Float64(x) for x ∈ fld.xrange])
    append!(rknots, xs[2:end-1])
  end
  sort!(unique!(rknots))
end

# apply the transmission correction to the spreading rate p when a ray crosses a
# soundspeed-gradient discontinuity at depth z_k: integrating the delta in
# dp/ds = -∂²c ξ² q through dz/ds = c ζ gives Δp = -q ξ² Δ(∂c/∂z) / (c |ζ|),
# independent of travel direction. Diverges at grazing (ray turning at the knot);
# skipped for |sin θ| ≤ 1e-3, consistent with _reflect_pq. If the crossing point
# lies at or beyond a boundary (a knot depth can coincide with the seabed or
# surface, e.g. at a seamount top), the ray reflects rather than transmits:
# terminate with the corresponding boundary event instead, as the restart could
# otherwise land outside the water column and tunnel through the boundary.
function _knot_jump!(i, (z_k, g₋, g₊, c_k), a, b, evref)
  u = i.u   # (r, z, ξ, ζ, t, p, q)
  if u[2] ≥ value(a, (u[1], 0.0, 0.0))
    evref[] = 1
    terminate!(i)
    return
  elseif u[2] ≤ -value(b, (u[1], 0.0, 0.0))
    evref[] = 2
    terminate!(i)
    return
  end
  ξ, ζ, q = u[3], u[4], u[7]
  abs(ζ) * c_k ≤ 1e-3 && return
  Δp = -q * ξ^2 * (g₊ - g₋) / (c_k * abs(ζ))
  i.u = Base.setindex(u, u[6] + Δp, 6)
  nothing
end

# prepare to trace a ray segment
function _prepare_trace(T, pm)
  ss = pm.env.soundspeed
  if ss isa SampledFieldZ && ss.interp isa UnderwaterAcoustics.Linear
    # where the water column does not extend beyond the sampled domain, evaluate
    # strictly inside the domain and continue the end segment linearly past its
    # edge: the field's flat extrapolation zeroes ∂c at/outside the outermost
    # knots, and when the bathymetry coincides with the last knot this biases
    # the curvature of every bottom-reflected ray (issue #29); clamping also
    # avoids Dual-indexed evaluation exactly at the domain edge
    # (Interpolations.jl #645). Where the water column extends beyond the
    # domain, the flat-extrapolated region is genuinely reachable and is left
    # untouched (the edge knot is then handled as a gradient discontinuity).
    z1, z2 = Float64.(extrema(ss.zrange))
    zlo = z1 ≤ -maximum(pm.env.bathymetry) + 1e-3 ? z1 + 1e-6 : -Inf
    zhi = z2 ≥ minimum(pm.env.altimetry) - 1e-3 ? z2 - 1e-6 : Inf
    c = z -> begin
      zc = clamp(z, zlo, zhi)
      value(ss, zc) + ∂(ss, zc, :z) * (z - zc)
    end
    ∂c = z -> ∂(ss, clamp(z, zlo, zhi), :z)
    ∂²c = z -> ∂(ss, clamp(z, zlo, zhi), :z, :z)
  else
    c = z -> value(ss, z)
    ∂c = z -> ∂(ss, z, :z)
    ∂²c = z -> ∂(ss, z, :z, :z)
  end
  ODEProblem{false}(_∂u, zeros(SVector{7,T}), (zero(T), one(T)), (c, ∂c, ∂²c))
end

# trace a ray starting at (r0, z0) with angle θ until it hits the surface,
# bottom, a scatterer, or reaches the domain limits. T is the type to use for
# computations, and p0, q0 are the initial spreading parameters (default to 1/c₀
# and 0, respectively). ds controls the spacing of points along the ray
# (0 = only events). guard suppresses scatterer detection for a short arc length
# after a scatterer bounce; crossbuf collects receiver-range crossings (event 6)
# and evref reports which terminal event fired (0 = none).
# Returns a 2-tuple of vectors containing distances along the ray, and
# (r, z, ξ, ζ, t, p, q) values at each point.

# a segment starting exactly on a knot makes the knot-crossing event fire at
# s = 0 (its root coincides with the initial condition), which collides with
# the knot step-size limiter and aborts the solve with dt → 0; start a hair
# off the knot, in the ray's initial vertical direction (Dual-safe: the
# comparisons use primal values, and the ternary keeps z0's own type)
_nudge_off_knots(z0, θ, ::Nothing, εz) = z0
function _nudge_off_knots(z0, θ, knots, εz)
  isempty(knots.list) && return z0
  minimum(k -> abs(z0 - k[1]), knots.list) < εz || return z0
  z0 + (sin(θ) < 0 ? -εz : εz) * one(z0)
end

function _trace_segment(T, prob, pm, r0, z0, θ, rmax, ds, p0, q0, guard, crossbuf, evref, r_rx, hmax, knots, rknots)
  a = pm.env.altimetry
  b = pm.env.bathymetry
  ss = pm.scatterers
  bs = pm.backscatter
  z0 = _nudge_off_knots(z0, θ, knots, 1e-3 * pm.atol)
  cᵥ = prob.p[1](z0)
  u0 = SA[convert(T, r0), z0, cos(θ)/cᵥ, sin(θ)/cᵥ, zero(T), p0, q0]
  # legacy tspan formula assumes forward-going rays (|θ| < π/2); with backscatter
  # or scatterers, use a geometry-based bound instead (segments end at events)
  rmin = one(T) * pm.rmin
  span = bs || ss !== nothing ? one(T) * pm.rugosity * (rmax - rmin + 4hmax) : one(T) * pm.rugosity * (rmax-r0)/cos(θ)
  tspan = (zero(T), span)
  prob = remake(prob; u0, tspan)
  δ = one(T) * pm.atol
  nk = knots === nothing ? 0 : length(knots.list)
  cb = VectorContinuousCallback(
    (out, u, s, i) -> _check_ray!(out, u, s, i, a, b, rmax, rmin, bs, ss, δ, guard, r_rx, knots, rknots),
    (i, ndx) -> _ray_event!(i, ndx, crossbuf, evref, knots, a, b), 6 + nk + length(rknots);
    rootfind = true
  )
  if ss === nothing
    cbs = cb
  else
    # anti-tunneling: since the independent variable is arc length and distance
    # is 1-Lipschitz, dt ≤ 0.9 × distance guarantees a step cannot jump over a
    # scatterer without the event function bracketing the surface (sphere tracing)
    limiter = StepsizeLimiter((u, pp, s) -> max(_mindist(ss, u[1], u[2]), δ); safety_factor=9//10, cached_dtcache=zero(T))
    cbs = CallbackSet(limiter, cb)
  end
  if knots !== nothing
    # the ODE error estimator is blind to the ∂c discontinuity at SSP knots, so
    # an integration step straddling a knot carries an unestimated error; shrink
    # steps approaching a knot so that the straddling step is at most δ long.
    # Over an arc s, |Δz| ≤ s (sinθ₀ + κ s) with κ = max |dθ/ds| = max |∂c|/c,
    # so a step of the positive root of s (sinθ₀ + κ s) = d cannot cross a knot
    # at vertical distance d (sinθ₀ overestimated via cmax for safety)
    klist, κ, cmax = knots.list, knots.κ, knots.cmax
    klimiter = StepsizeLimiter(
      (u, pp, s) -> begin
        d = minimum(k -> abs(u[2] - k[1]), klist)
        sinθ = abs(u[4]) * cmax
        max((√(sinθ^2 + 4κ * d) - sinθ) / (2κ), δ)
      end; safety_factor=9//10, cached_dtcache=zero(T))
    cbs = CallbackSet(klimiter, cbs)
  end
  # reltol must track solver_tol: the ODE-solver default (1e-3) would dominate
  # abstol for states of O(10+) (depths, ranges) and leave a mm-scale geometry
  # error floor that no abstol can remove
  if ds ≤ 0
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, reltol=pm.solver_tol, save_everystep=false, callback=cbs)
  else
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, reltol=pm.solver_tol, saveat=ds, callback=cbs)
  end
  (soln.t, soln.u)
end

# dynamic ray tracing correction to the spreading rate p (= qp) for reflection
# off a (possibly curved) boundary in a medium with a soundspeed gradient
# [COA (3.123-3.128); Müller, Geoph. J. R.A.S. 79 (1984); BELLHOP Reflect2D].
# t̂ is the incident unit ray tangent, n̂ the boundary unit normal pointing into
# the water, κ the boundary curvature (positive where the boundary is convex
# toward the water), and gc the soundspeed gradient (∂c/∂x, ∂c/∂z) at the
# reflection point. The correction diverges at grazing incidence and is
# skipped for sinθg ≤ 1e-3, consistent with the eigenray search tolerances.
function _reflect_pq(qp, q, c, t̂, n̂, κ, gc)
  Th = dot(t̂, n̂) / c                    # normal component of the 1/c-scaled tangent
  abs(Th) * c ≤ 1e-3 && return qp
  Tg = (t̂[1] * n̂[2] - t̂[2] * n̂[1]) / c  # tangential component (t̂b = (n̂z, -n̂x))
  t̂′ = t̂ - 2 * dot(t̂, n̂) * n̂
  rayn = SA[-t̂[2], t̂[1]]                # ray-centred unit normals
  rayn′ = SA[t̂′[2], -t̂′[1]]             # (reflected frame has opposite sense)
  cnjump = -dot(gc, rayn′ - rayn)
  csjump = -dot(gc, t̂′ - t̂)
  RM = Tg / Th
  RN = -2κ / (c^2 * Th) + RM * (2 * cnjump - RM * csjump) / c^2
  qp + q * RN
end

# trace a ray starting at tx1 with angle θ until it reaches rmax. ds controls the
# spacing of points along the ray. If paths=true, also returns the ray path as a
# vector of (x, y, z) points. If aux=true, also returns auxiliary info as a vector
# of (s, q, t, A) tuples at each point. If r_rx is given (not NaN), receiver-range
# crossings are recorded and returned as the third element of the result.
function _trace(pm::RaySolver, tx1::AbstractAcousticSource, θ, rmax, ds=0.0; paths=true, aux=false, r_rx=NaN)
  θ₀ = θ
  f = frequency(tx1)
  mid_temp = (minimum(pm.env.temperature) + maximum(pm.env.temperature)) / 2
  hmax = maximum(pm.env.bathymetry)
  zmin = -hmax
  ϵ = one(zmin) * pm.solver_tol
  p = location(tx1)
  T = promote_type(env_type(pm.env), eltype(p), typeof(f), typeof(θ), typeof(rmax), _scat_eltype(pm.scatterers))
  raypath = Array{@NamedTuple{x::T,y::T,z::T}}(undef, 0)
  aux_info = Array{Tuple{T,T,T,Complex{T}}}(undef, 0)
  crossings = Array{@NamedTuple{z::T,t::T,θ::T,q::T,A::Complex{T},dir::Int,D::T,ns::Int,nb::Int,nsc::Int,pathlen::Int}}(undef, 0)
  crossbuf = Array{Tuple{T,SVector{7,T}}}(undef, 0)
  evref = Ref(0)
  rxr = convert(T, r_rx)
  prob = _prepare_trace(T, pm)
  cf, ∂cf, _ = prob.p    # domain-safe soundspeed and gradient (see _prepare_trace)
  c₀ = cf(p.z)
  A = one(Complex{T})   # phasor
  s = 0                 # surface bounces
  b = 0                 # bottom bounces
  sc = 0                # scatterer bounces
  t = zero(T)           # time along ray
  D = zero(T)           # distance along ray
  q = zero(T)           # spreading factor
  qp = one(T) / c₀      # spreading rate
  guard = zero(T)       # arc length over which scatterer detection is suppressed
  knots = _gradient_discontinuities(pm.env)
  rknots = _boundary_range_knots(pm.env)
  while true
    empty!(crossbuf)
    evref[] = 0
    npath0 = length(raypath)
    svec, u = _trace_segment(T, prob, pm, p[1], p[3], θ, rmax, ds, qp, q, guard, crossbuf, evref, rxr, hmax, knots, rknots)
    guard = zero(T)
    r, z, ξ, ζ, dt, qp, q = u[end]
    dD = svec[end]
    θ = atan(ζ, ξ)
    oq = u[1][7]
    kmah = 0
    for i ∈ eachindex(u)
      sign(oq) * sign(u[i][7]) == -1 && (kmah += 1)
      u[i][7] != 0.0 && (oq = u[i][7])
      paths && push!(raypath, xyz(u[i][1], 0.0, u[i][2]))
    end
    if aux
      for i ∈ eachindex(u)
        push!(aux_info, (D + svec[i], u[i][7], t + u[i][5], A))
      end
    end
    for (sx, ux) ∈ crossbuf     # receiver-range crossings within this segment
      kx = 0                    # KMAH count up to the crossing point
      npts = 0
      oqx = u[1][7]
      for i ∈ eachindex(u)
        svec[i] > sx && break
        npts += 1
        sign(oqx) * sign(u[i][7]) == -1 && (kx += 1)
        u[i][7] != 0.0 && (oqx = u[i][7])
      end
      push!(crossings, (; z=ux[2], t=t+ux[5], θ=atan(ux[4], ux[3]), q=ux[7],
        A=A*cis(-π/2*kx), dir=Int(sign(ux[3])), D=D+sx, ns=s, nb=b, nsc=sc,
        pathlen=paths ? npath0+npts : 0))
    end
    t += dt
    D += dD
    A *= cis(-π/2 * kmah)                 # [COA §3.4.1 (KMAH correction)]
    zmin = -value(pm.env.bathymetry, (r, 0.0, 0.0))
    zmax = value(pm.env.altimetry, (r, 0.0, 0.0))
    p = (r, 0.0, clamp(z, zmin+ϵ, zmax-ϵ))
    ev = evref[]
    if r ≥ rmax - 1e-3
      break
    elseif ev == 5                        # hit a scatterer
      sc += 1
      sca = pm.scatterers.items[_nearest(pm.scatterers, r, z)]
      hit = UnderwaterAcoustics.boundary_projection(sca.shape, (x=r, y=zero(r), z=z))
      n̂ = SA[hit.normal.x, hit.normal.z]
      t̂ = SA[cos(θ), sin(θ)]
      dp = dot(t̂, n̂)
      t̂′ = t̂ - 2dp * n̂                    # specular reflection off local tangent plane
      θ = atan(t̂′[2], t̂′[1])
      sinθg = clamp(abs(dp), zero(dp), one(dp))
      ρ = value(pm.env.density, (r, 0.0, z))
      c = cf(z)
      A *= reflection_coef(sca.boundary, f, asin(sinθg), ρ, c)
      gz = ∂cf(z)
      qp = _reflect_pq(qp, q, c, t̂, n̂, hit.curvature, SA[zero(gz), gz])
      guard = 2 * one(T) * pm.atol
    elseif ev == 1 && isapprox(z, zmax; atol=1e-3)   # hit the surface
      s += 1
      ρ = value(pm.env.density, (r, 0.0, 0.0))
      c = cf(zero(z))
      A *= reflection_coef(pm.env.surface, f, π/2 - θ, ρ, c)
      a′ = ∂(pm.env.altimetry, (r, 0.0, 0.0), :x)
      a″ = ∂(pm.env.altimetry, (r, 0.0, 0.0), :x, :x)
      m = √(1 + a′^2)
      κ = a″ / m^3                       # positive: boundary convex into the water
      n̂ = SA[a′, -one(a′)] / m           # unit normal pointing into the water
      gz = ∂cf(z)
      qp = _reflect_pq(qp, q, c, SA[cos(θ), sin(θ)], n̂, κ, SA[zero(gz), gz])
      θ = -θ + 2atan(a′)
    elseif ev == 2 && isapprox(z, zmin; atol=1e-3)   # hit the bottom
      b += 1
      ρ = value(pm.env.density, (r, 0.0, z))
      c = cf(z)
      A *= reflection_coef(pm.env.seabed, f, π/2 + θ, ρ, c)
      d′ = ∂(pm.env.bathymetry, (r, 0.0, 0.0), :x)
      d″ = ∂(pm.env.bathymetry, (r, 0.0, 0.0), :x, :x)
      m = √(1 + d′^2)
      κ = d″ / m^3                       # positive: boundary convex into the water
      n̂ = SA[d′, one(d′)] / m            # unit normal pointing into the water
      gz = ∂cf(z)
      qp = _reflect_pq(qp, q, c, SA[cos(θ), sin(θ)], n̂, κ, SA[zero(gz), gz])
      θ = -θ - 2atan(d′)
    else
      break
    end
    if pm.backscatter
      abs(A) * absorption(f, D, pm.env.salinity, mid_temp, hmax/2) / abs(q) < pm.min_amplitude && break
    else
      abs(θ) ≥ π/2 && break
      abs(A)/abs(q) < pm.min_amplitude && break
    end
  end
  cₛ = cf(p[3])
  A *= √abs(cₛ * cos(θ₀) / (p[1] * c₀ * q))                   # [COA (3.65)]
  A *= absorption(f, D, pm.env.salinity, mid_temp, -zmin/2)   # nominal absorption
  RayArrival(t, A, s, b, θ₀, -θ, raypath), aux_info, crossings
end

### eigenray search (signature-grouped; used by both forward and backscatter paths)

# crossings are grouped by "signature" — direction of travel + bounce history —
# so that the depth error Δz(θ) restricted to one signature is continuous in θ
# and the standard bracketing/root-finding machinery applies
_sigbase(cr) = (cr.dir, cr.ns, cr.nb, cr.nsc)

# signature of each crossing: (dir, ns, nb, nsc, ordinal within that history)
function _sigs(crossings)
  counts = Dict{NTuple{4,Int},Int}()
  map(crossings) do cr
    base = _sigbase(cr)
    k = get(counts, base, 0) + 1
    counts[base] = k
    (base..., k)
  end
end

# index of the crossing matching a signature, or nothing
function _findsig(crossings, sig)
  base = sig[1:4]
  k = 0
  for i ∈ eachindex(crossings)
    if _sigbase(crossings[i]) == base
      k += 1
      k == sig[5] && return i
    end
  end
  nothing
end

# depth error at the receiver range for the crossing matching sig; NaN when the
# launch angle is outside the fan (a diverging Newton iteration can probe
# arbitrary angles, and without backscatter a beyond-vertical angle would give a
# negative integration span)
function _Δz_sig(ϕ, (pm, tx1, rdom, r_rx, z, sig))
  ϕv = ForwardDiff.value(ϕ)
  if !(pm.backscatter ? -π ≤ ϕv ≤ π : abs(ϕv) < π/2)
    return convert(promote_type(typeof(ϕ), typeof(z)), NaN)
  end
  crs = _trace(pm, tx1, ϕ, rdom; paths=false, r_rx)[3]
  i = _findsig(crs, sig)
  i === nothing ? convert(promote_type(typeof(ϕ), typeof(z)), NaN) : crs[i].z - z
end

# trace a ray and construct a RayArrival from the crossing matching sig,
# applying the COA (3.65) amplitude formula and nominal absorption at the
# crossing point; nothing if the signature is absent for this launch angle
function _crossing_arrival(pm, tx1, θ, rdom, r_rx, sig, ds; paths=true)
  f = frequency(tx1)
  mid_temp = (minimum(pm.env.temperature) + maximum(pm.env.temperature)) / 2
  hmax = maximum(pm.env.bathymetry)
  c₀ = value(pm.env.soundspeed, location(tx1))
  arr, _, crs = _trace(pm, tx1, θ, rdom, ds; paths, r_rx)
  i = _findsig(crs, sig)
  i === nothing && return nothing
  cr = crs[i]
  cₛ = value(pm.env.soundspeed, (r_rx, 0.0, cr.z))
  A = cr.A * √abs(cₛ * cos(θ) / (r_rx * c₀ * cr.q))           # [COA (3.65)]
  A *= absorption(f, cr.D, pm.env.salinity, mid_temp, hmax/2) # nominal absorption
  path = paths ? [arr.path[1:cr.pathlen]; xyz(oftype(cr.z, r_rx), zero(cr.z), cr.z)] : arr.path[1:0]
  RayArrival(cr.t, A, cr.ns, cr.nb, θ, -cr.θ, path)
end

function _arrivals_backscatter(pm::RaySolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true, ztol=0.1)
  p2 = location(rx)
  r_rx = p2.x
  hmax = maximum(pm.env.bathymetry)
  rdom = pm.rmax > 0 ? pm.rmax : r_rx + 0.1
  nbeams = pm.nbeams
  if nbeams == 0
    R = max(abs(r_rx - location(tx).x), hmax)
    nbeams = ceil(Int, 16 * (pm.max_angle - pm.min_angle) / atan(hmax, R))
  end
  ds = pm.ds ≤ 0 ? minimum(pm.env.bathymetry) / 10 : pm.ds
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  T1 = promote_type(env_type(pm.env), eltype(location(tx)), typeof(p2.x), typeof(p2.z))
  T2 = typeof(RayArrival(zero(T1), zero(Complex{T1}), 0, 0, zero(T1), zero(T1),
    Array{@NamedTuple{x::T1,y::T1,z::T1}}(undef, 0)))
  erays = _fan_search(T2, pm, tx, collect(T1, θ), T1(rdom), T1(r_rx), T1(p2.z), ds, ztol, paths)
  if pm.rmin < 0    # also insonify the left halfspace with a mirrored fan
    append!(erays, _fan_search(T2, pm, tx, collect(T1, rem2pi.(π .- θ, RoundNearest)), T1(rdom), T1(r_rx), T1(p2.z), ds, ztol, paths))
  end
  sort!(erays; by=Base.Fix2(getfield, :t))
end

# eigenray search over one launch-angle fan: find receiver-range crossings for
# each launch angle, group by signature, and root-find on the per-signature
# depth error (exact roots, brackets, and near-turning-point double roots).
# With edges=true, signature support edges are additionally refined by
# bisection (issue #31: recovers families whose roots hide between the last
# supported fan sample and the tangency angle behind a blocking feature); off
# for backscatter, where the near-tangent grazing arrivals it recovers make
# the echo field noisy w.r.t. scatterer geometry.
function _fan_search(::Type{T2}, pm::RaySolver, tx::AbstractAcousticSource, θ, rdom, r_rx, zrx, ds, ztol, paths; edges=false) where T2
  T1 = eltype(θ)
  # the interval root-solve tolerance is in launch-angle space, but eigenrays are
  # accepted on the depth residual (ztol, in m); since dΔz/dθ grows ∝ range, the
  # angle tolerance must shrink ∝ 1/range or all roots at long range carry
  # residuals ≫ ztol and are rejected
  θtol = pm.atol / max(ForwardDiff.value(abs(r_rx)), 1.0)   # plain float: solver control, not differentiated
  ntasks = _ntasks(pm)
  crs_all = _tmap(θ1 -> _trace(pm, tx, θ1, rdom; paths=false, r_rx)[3], θ, ntasks)
  allsigs = Set{NTuple{5,Int}}()
  foreach(crs -> union!(allsigs, _sigs(crs)), crs_all)
  erays = T2[]
  for sig ∈ allsigs
    err = map(crs_all) do crs
      i = _findsig(crs, sig)
      i === nothing ? T1(NaN) : T1(crs[i].z - zrx)
    end
    prob1 = IntervalNonlinearProblem{false}(_Δz_sig, T1.((θ[1], θ[1])), (pm, tx, rdom, r_rx, zrx, sig))
    prob2 = NonlinearProblem{false}(_Δz_sig, T1(θ[1]), (pm, tx, rdom, r_rx, zrx, sig))
    results = Vector{T2}[T2[] for _ ∈ 1:length(θ)]
    _tforeach(length(θ), ntasks) do i
      if isapprox(err[i], 0; atol=pm.atol)
        # crossing already at receiver
        v = _crossing_arrival(pm, tx, θ[i], rdom, r_rx, sig, ds; paths)
        v === nothing || push!(results[i], v)
      elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) < 0
        # crossings bracket the receiver, so find a root in between...
        soln = solve(remake(prob1; tspan=T1.(_ordered(θ[i-1], θ[i]))); abstol=θtol)
        if successful_retcode(soln.retcode) && abs(soln.resid) < ztol
          v = _crossing_arrival(pm, tx, soln.u, rdom, r_rx, sig, ds; paths)
          v === nothing || push!(results[i], v)
        end
      elseif i > 2 && _isnearzero(err[i-2], err[i-1], err[i])
        # at a turning point, so potentially two roots between i-2 and i, try and find both...
        lims = _ordered(θ[i-2], θ[i])
        soln = solve(remake(prob2; u0=T1(θ[i-1])); abstol=pm.atol)
        if successful_retcode(soln.retcode) && abs(soln.resid) < ztol && lims[1] < soln.u < lims[2]
          θ₁ = soln.u
          v = _crossing_arrival(pm, tx, θ₁, rdom, r_rx, sig, ds; paths)
          v === nothing || push!(results[i], v)
          soln = solve(remake(prob2; u0=T1(θ[i-1] < θ₁ ? lims[2] : lims[1])); abstol=pm.atol)
          if successful_retcode(soln.retcode) && abs(soln.resid) < ztol &&
             (θ[i-1] < θ₁ ? θ₁ < soln.u < lims[2] : lims[1] < soln.u < θ₁)
            v = _crossing_arrival(pm, tx, soln.u, rdom, r_rx, sig, ds; paths)
            v === nothing || push!(results[i], v)
          end
        end
      end
      # support-edge refinement: a signature's θ-support ends where the ray
      # family becomes tangent to a blocking feature; a root can hide between
      # the last supported fan sample and that tangency angle (with no second
      # sample to bracket it). Bisect toward the edge, and root-find on any
      # sign change encountered on the way.
      if edges && !isnan(err[i]) && !isapprox(err[i], 0; atol=pm.atol)
        for j ∈ (i-1, i+1)
          1 ≤ j ≤ length(θ) && isnan(err[j]) || continue
          θr = _edge_root(θ[i], θ[j], err[i], θtol, (pm, tx, rdom, r_rx, zrx, sig), ztol)
          if θr !== nothing
            v = _crossing_arrival(pm, tx, θr, rdom, r_rx, sig, ds; paths)
            v === nothing || push!(results[i], v)
          end
        end
      end
    end
    foreach(r -> append!(erays, r), results)
  end
  erays
end

# bisect from a supported launch angle θa toward an unsupported one θb (where
# the signature is absent) and return a root of the per-signature depth error
# if a sign change is found before the support edge; nothing otherwise. The
# θ-support of a signature can be porous near tangency (a grazing ray's
# hit/miss of the blocking feature flips with integration noise), which
# defeats the standard bracketing solvers, so the search runs a NaN-aware
# bisection on primal values; for dual-numbered inputs the root's derivatives
# are then restored via the implicit function theorem (matching what
# NonlinearSolve does for a clean bracket).
function _edge_root(θa, θb, fa, θtol, params, ztol)
  pf = ϕ -> ForwardDiff.value(_Δz_sig(ϕ, params))
  lo, hi = ForwardDiff.value(θa), ForwardDiff.value(θb)
  flo = ForwardDiff.value(fa)
  found = false
  local θ2, f2
  while abs(hi - lo) > θtol
    m = (lo + hi) / 2
    fm = pf(m)
    if isnan(fm)
      hi = m
    elseif sign(fm) * sign(flo) < 0
      θ2, f2 = m, fm
      found = true
      break
    else
      lo = m
    end
  end
  found || return nothing
  θ1, f1 = lo, flo
  # NaN-aware bisection on [θ1, θ2]; NaN midpoints are replaced by a nearby
  # supported angle (probes at ±j/20 of the interval), so a porous support
  # only slows the shrink instead of killing it
  while abs(θ2 - θ1) > θtol
    m = (θ1 + θ2) / 2
    fm = pf(m)
    if isnan(fm)
      ok = false
      for j ∈ 1:9, s ∈ (1.0, -1.0)
        x = m + s * j / 20 * (θ2 - θ1)
        fx = pf(x)
        isnan(fx) && continue
        m, fm = x, fx
        ok = true
        break
      end
      ok || break
    end
    if sign(fm) == sign(f1)
      θ1, f1 = m, fm
    else
      θ2, f2 = m, fm
    end
  end
  θr, fr = abs(f1) < abs(f2) ? (θ1, f1) : (θ2, f2)
  abs(fr) < ztol || return nothing
  T1 = typeof(θa)
  T1 <: ForwardDiff.Dual || return θr
  # implicit function theorem: dθ*/dp = -(∂f/∂p) / (∂f/∂θ) at the root
  slope = (f2 - f1) / (θ2 - θ1)
  fD = _Δz_sig(convert(T1, θr), params)
  isnan(ForwardDiff.value(fD)) && return convert(T1, θr)
  convert(T1, θr) - (fD - ForwardDiff.value(fD)) / slope
end
