import LinearAlgebra: norm, dot
import OrdinaryDiffEq: ODEProblem, VectorContinuousCallback, CallbackSet, Tsit5, solve
import OrdinaryDiffEqRosenbrock: Rosenbrock23
import DiffEqCallbacks: StepsizeLimiter
import NonlinearSolve: IntervalNonlinearProblem, NonlinearProblem
import SciMLBase: successful_retcode, terminate!, remake
import ForwardDiff: derivative
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
- `min_amplitude`: minimum ray amplitude to track (default: 1e-5)
- `solver`: differential equation solver (default: nothing, auto)
- `solver_tol`: solver tolerance (default: 1e-4)
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
  nbeams = pm.nbeams
  if nbeams == 0
    p1 = location(tx)
    R = abs(p2.x - p1.x)
    h = maximum(pm.env.bathymetry)
    nbeams = ceil(Int, 16 * (pm.max_angle - pm.min_angle) / atan(h, R))
  end
  ds = pm.ds ≤ 0 ? minimum(pm.env.bathymetry) / 10 : pm.ds
  ntasks = _ntasks(pm)
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  rays = _tmap(θ1 -> _trace(pm, tx, θ1, p2.x)[1], θ, ntasks)
  err = map(rays) do ray
    p3 = ray.path[end]
    isapprox(p3.x, p2.x; atol=pm.atol) && isapprox(p3.y, p2.y; atol=pm.atol) ? p3.z - p2.z : NaN
  end
  T1 = promote_type(env_type(pm.env), eltype(location(tx)), typeof(p2.x), typeof(p2.z))
  T2 = T1 == eltype(θ) ? eltype(rays) : typeof(_trace(pm, tx, T1(θ[1]), p2.x; paths)[1])
  prob1 = IntervalNonlinearProblem{false}(_Δz, T1.((θ[1], θ[1])), (pm, tx, p2.x, p2.z))
  prob2 = NonlinearProblem{false}(_Δz, T1(θ[1]), (pm, tx, p2.x, p2.z))
  erays = [T2[] for _ ∈ 1:length(θ)]
  _tforeach(length(θ), ntasks) do i
    if isapprox(err[i], 0; atol=pm.atol)
      # ray already at receiver
      v = T1 == eltype(θ) ? rays[i] : _trace(pm, tx, T1(θ[i]), p2.x; paths)[1]
      push!(erays[i], v)
    elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) < 0
      # rays bracket the receiver, so find a root in between...
      soln = solve(remake(prob1; tspan=T1.(_ordered(θ[i-1], θ[i]))); abstol=pm.atol)
      if successful_retcode(soln.retcode) && abs(soln.resid) < ztol
        push!(erays[i], _trace(pm, tx, soln.u, p2.x, ds; paths)[1])
      end
    elseif i > 2 && _isnearzero(err[i-2], err[i-1], err[i])
      # at a turning point, so potentially two roots between i-2 and i, try and find both...
      lims = _ordered(θ[i-2], θ[i])
      soln = solve(remake(prob2; u0=T1(θ[i-1])); abstol=pm.atol)
      if successful_retcode(soln.retcode) && abs(soln.resid) < ztol && lims[1] < soln.u < lims[2]
        θ₁ = soln.u
        push!(erays[i], _trace(pm, tx, θ₁, p2.x, ds; paths)[1])
        soln = solve(remake(prob2; u0=T1(θ[i-1] < θ₁ ? lims[2] : lims[1])); abstol=pm.atol)
        if successful_retcode(soln.retcode) && abs(soln.resid) < ztol &&
           (θ[i-1] < θ₁ ? θ₁ < soln.u < lims[2] : lims[1] < soln.u < θ₁)
          push!(erays[i], _trace(pm, tx, soln.u, p2.x, ds; paths)[1])
        end
      end
    end
  end
  sort!(reduce(vcat, erays); by=Base.Fix2(getfield, :t))
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
    nbeams = clamp(ceil(Int, 8 * (pm.max_angle - pm.min_angle) * R / Δz), 100, 1000)
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
        4 * max(abs(q1), abs(q2)) * δθ                  # size of neighborhood
      )
      for j ∈ xi, k ∈ zi                                # loop over the subset
        rx = rxs[j,k].pos
        rxpos = SA[rx.x, rx.z]
        α = dot(rxpos - pos1, vec12) / vec12_mag2
        0 ≤ α < 1 || continue                           # rx outside of ray segment
        q = q1 + α * (q2 - q1)                          # spreading factor at cpa
        W = abs(q * δθ)                                 # beam width at cpa [COA (3.74)]
        cpa = pos1 + α * vec12                          # closest point of approach
        cpa[1] > 0 || continue                          # no deposit at/behind the r=0 axis
        n = norm(cpa - rxpos)                           # normal distance from ray to rx
        n < 4W || continue                              # rx too far from ray segment
        A = C1 * sqrt(C2 / (cpa[1] * W))                # [COA (3.76)]
        if mode === :coherent
          t = t1 + α * (t2 - t1)                        # time at cpa
          P = A * exp(-(n / W)^2) * cispi(2f * t)       # [COA (3.72), COA (3.75)]
        else
          P = complex(abs2(A * exp(-(n / W)^2)), 0.0)
        end
        afld[j,k] += P
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
# crossing (recorded, not terminal). Inert slots hold one(u[1]) so they never fire
# and preserve the element type (dual numbers).
function _check_ray!(out, u, s, integrator, a, b, rmax, rmin, backscatter, ss, δ, guard, r_rx)
  out[1] = value(a, (u[1], 0.0, 0.0)) - u[2]    # surface reflection
  out[2] = u[2] + value(b, (u[1], 0.0, 0.0))    # bottom reflection
  out[3] = rmax - u[1]                          # maximum range
  out[4] = backscatter ? u[1] - rmin + δ : u[3] # exit domain left / ray turned back
  out[5] = ss === nothing || s < guard ? one(u[1]) : _mindist(ss, u[1], u[2]) - δ
  out[6] = isnan(r_rx) ? one(u[1]) : u[1] - r_rx
end

# handle a fired ray event: event 6 (receiver-range crossing) is recorded and
# integration continues; any other event terminates the segment and its index is
# reported via evref. Depending on the DiffEqBase version, ndx is either the
# event index or a vector of simultaneous event flags (0 = not fired, ±1 = fired)
function _ray_event!(i, ndx::Integer, crossbuf, evref)
  if ndx == 6
    push!(crossbuf, (i.t, i.u))
  else
    evref[] = Int(ndx)
    terminate!(i)
  end
end

function _ray_event!(i, ndx::AbstractVector, crossbuf, evref)
  length(ndx) ≥ 6 && ndx[6] != 0 && push!(crossbuf, (i.t, i.u))
  for k ∈ 1:min(length(ndx), 5)
    if ndx[k] != 0
      evref[] = k
      terminate!(i)
      break
    end
  end
end

# prepare to trace a ray segment
function _prepare_trace(T, pm)
  c = z -> value(pm.env.soundspeed, z)
  ∂c = z -> derivative(c, z)
  ∂²c = z -> derivative(∂c, z)
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
function _trace_segment(T, prob, pm, r0, z0, θ, rmax, ds, p0, q0, guard, crossbuf, evref, r_rx, hmax)
  a = pm.env.altimetry
  b = pm.env.bathymetry
  ss = pm.scatterers
  bs = pm.backscatter
  cᵥ = value(pm.env.soundspeed, z0)
  u0 = SA[convert(T, r0), z0, cos(θ)/cᵥ, sin(θ)/cᵥ, zero(T), p0, q0]
  # legacy tspan formula assumes forward-going rays (|θ| < π/2); with backscatter
  # or scatterers, use a geometry-based bound instead (segments end at events)
  rmin = one(T) * pm.rmin
  span = bs || ss !== nothing ? one(T) * pm.rugosity * (rmax - rmin + 4hmax) : one(T) * pm.rugosity * (rmax-r0)/cos(θ)
  tspan = (zero(T), span)
  prob = remake(prob; u0, tspan)
  δ = one(T) * pm.atol
  cb = VectorContinuousCallback(
    (out, u, s, i) -> _check_ray!(out, u, s, i, a, b, rmax, rmin, bs, ss, δ, guard, r_rx),
    (i, ndx) -> _ray_event!(i, ndx, crossbuf, evref), 6;
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
  if ds ≤ 0
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, save_everystep=false, callback=cbs)
  else
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, saveat=ds, callback=cbs)
  end
  (soln.t, soln.u)
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
  c₀ = value(pm.env.soundspeed, p)
  T = promote_type(env_type(pm.env), eltype(p), typeof(f), typeof(θ), typeof(rmax), _scat_eltype(pm.scatterers))
  raypath = Array{@NamedTuple{x::T,y::T,z::T}}(undef, 0)
  aux_info = Array{Tuple{T,T,T,Complex{T}}}(undef, 0)
  crossings = Array{@NamedTuple{z::T,t::T,θ::T,q::T,A::Complex{T},dir::Int,D::T,ns::Int,nb::Int,nsc::Int,pathlen::Int}}(undef, 0)
  crossbuf = Array{Tuple{T,SVector{7,T}}}(undef, 0)
  evref = Ref(0)
  rxr = convert(T, r_rx)
  prob = _prepare_trace(T, pm)
  A = one(Complex{T})   # phasor
  s = 0                 # surface bounces
  b = 0                 # bottom bounces
  sc = 0                # scatterer bounces
  t = zero(T)           # time along ray
  D = zero(T)           # distance along ray
  q = zero(T)           # spreading factor
  qp = one(T) / c₀      # spreading rate
  guard = zero(T)       # arc length over which scatterer detection is suppressed
  while true
    empty!(crossbuf)
    evref[] = 0
    npath0 = length(raypath)
    svec, u = _trace_segment(T, prob, pm, p[1], p[3], θ, rmax, ds, qp, q, guard, crossbuf, evref, rxr, hmax)
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
      c = value(pm.env.soundspeed, (r, 0.0, z))
      A *= reflection_coef(sca.boundary, f, asin(sinθg), ρ, c)
      # dynamic ray tracing correction for reflection off a curved boundary
      # [COA (3.123-3.128)]; skipped near grazing where the correction diverges
      sinθg > 1e-3 && (qp += q * 2 * hit.curvature / (c * sinθg))
      guard = 2 * one(T) * pm.atol
    elseif ev == 1 && isapprox(z, zmax; atol=1e-3)   # hit the surface
      s += 1
      ρ = value(pm.env.density, (r, 0.0, 0.0))
      c = value(pm.env.soundspeed, (r, 0.0, 0.0))
      A *= reflection_coef(pm.env.surface, f, π/2 - θ, ρ, c)
      α = atan(derivative(x -> value(pm.env.altimetry, (x, 0.0, 0.0)), r))
      θ = -θ + 2α
    elseif ev == 2 && isapprox(z, zmin; atol=1e-3)   # hit the bottom
      b += 1
      ρ = value(pm.env.density, (r, 0.0, z))
      c = value(pm.env.soundspeed, (r, 0.0, z))
      A *= reflection_coef(pm.env.seabed, f, π/2 + θ, ρ, c)
      α = atan(derivative(x -> value(pm.env.bathymetry, (x, 0.0, 0.0)), r))
      θ = -θ - 2α
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
  cₛ = value(pm.env.soundspeed, p)
  A *= √abs(cₛ * cos(θ₀) / (p[1] * c₀ * q))                   # [COA (3.65)]
  A *= absorption(f, D, pm.env.salinity, mid_temp, -zmin/2)   # nominal absorption
  RayArrival(t, A, s, b, θ₀, -θ, raypath), aux_info, crossings
end

_Δz(ϕ, (pm, tx1, rmax, z)) = _trace(pm, tx1, ϕ, rmax)[1].path[end].z - z

### backscatter eigenray search

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

# depth error at the receiver range for the crossing matching sig
function _Δz_sig(ϕ, (pm, tx1, rdom, r_rx, z, sig))
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
# depth error using the same 3-way logic as the forward eigenray search
function _fan_search(T2::Type, pm::RaySolver, tx::AbstractAcousticSource, θ, rdom, r_rx, zrx, ds, ztol, paths)
  T1 = eltype(θ)
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
    results = [T2[] for _ ∈ 1:length(θ)]
    _tforeach(length(θ), ntasks) do i
      if isapprox(err[i], 0; atol=pm.atol)
        # crossing already at receiver
        v = _crossing_arrival(pm, tx, θ[i], rdom, r_rx, sig, ds; paths)
        v === nothing || push!(results[i], v)
      elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) < 0
        # crossings bracket the receiver, so find a root in between...
        soln = solve(remake(prob1; tspan=T1.(_ordered(θ[i-1], θ[i]))); abstol=pm.atol)
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
    end
    foreach(r -> append!(erays, r), results)
  end
  erays
end
