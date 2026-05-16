import LinearAlgebra: norm, dot
import OrdinaryDiffEq: ODEProblem, VectorContinuousCallback, Tsit5, solve
import OrdinaryDiffEqRosenbrock: Rosenbrock23
import NonlinearSolve: IntervalNonlinearProblem, NonlinearProblem
import SciMLBase: successful_retcode, terminate!
import ForwardDiff: derivative
import StaticArrays: SA

Base.@kwdef struct RaySolver{T1,T2} <: AbstractRayPropagationModel
  env::T1
  nbeams::Int = 0
  min_angle::Float64 = -deg2rad(80)
  max_angle::Float64 = +deg2rad(80)
  ds::Float64 = 1.0
  atol::Float64 = 1e-4
  rugosity::Float64 = 1.5
  min_amplitude::Float64 = 1e-6
  solver::T2 = nothing
  solver_tol::Float64 = 1e-8
  function RaySolver(env, nbeams, min_angle, max_angle, ds, atol, rugosity, min_amplitude, solver, solver_tol)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ min_angle ≤ π/2 || error("min_angle should be between -π/2 and π/2")
    -π/2 ≤ max_angle ≤ π/2 || error("max_angle should be between -π/2 and π/2")
    min_angle < max_angle || error("max_angle should be more than min_angle")
    solver = something(solver, is_isovelocity(env) ? Tsit5() : Rosenbrock23())
    new{typeof(env),typeof(solver)}(env, nbeams, min_angle, max_angle, ds, atol, rugosity, min_amplitude, solver, solver_tol)
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
- `ds`: nominal spacing between ray points (default: 1 m)
- `atol`: absolute position tolerance (default: 0.0001 m)
- `rugosity`: TODO (default: 1.5)
- `min_amplitude`: minimum ray amplitude to track (default: 1e-5)
- `solver`: differential equation solver (default: nothing, auto)
- `solver_tol`: solver tolerance (default: 1e-4)
"""
RaySolver(env; kwargs...) = RaySolver(; env=env, kwargs...)

Base.show(io::IO, pm::RaySolver) = print(io, "RaySolver(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::RaySolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true, ztol=0.1)
  _check2d([tx], [rx])
  p2 = location(rx)
  nbeams = pm.nbeams
  if nbeams == 0
    p1 = location(tx)
    R = abs(p2.x - p1.x)
    h = maximum(pm.env.bathymetry)
    nbeams = ceil(Int, 16 * (pm.max_angle - pm.min_angle) / atan(h, R))
  end
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  rays = tmap(θ1 -> _trace(pm, tx, θ1, p2.x)[1], θ)
  err = map(rays) do ray
    p3 = ray.path[end]
    isapprox(p3.x, p2.x; atol=pm.atol) && isapprox(p3.y, p2.y; atol=pm.atol) ? p3.z - p2.z : NaN
  end
  T1 = promote_type(env_type(pm.env), eltype(location(tx)), typeof(p2.x), typeof(p2.z))
  T2 = T1 == eltype(θ) ? eltype(rays) : typeof(_trace(pm, tx, T1(θ[1]), p2.x; paths)[1])
  erays = Channel{T2}(length(θ))
  Threads.@threads for i ∈ 1:length(θ)
    if isapprox(err[i], 0; atol=pm.atol)
      # ray already at receiver
      v = T1 == eltype(θ) ? rays[i] : _trace(pm, tx, T1(θ[i]), p2.x; paths)[1]
      put!(erays, v)
    elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) < 0
      # rays bracket the receiver, so find a root in between...
      soln = solve(IntervalNonlinearProblem{false}(_Δz, T1.(_ordered(θ[i-1], θ[i])), (pm, tx, p2.x, p2.z)); abstol=pm.atol)
      if successful_retcode(soln.retcode) && abs(soln.resid) < ztol
        put!(erays, _trace(pm, tx, soln.u, p2.x, pm.ds; paths)[1])
      end
    elseif i > 2 && _isnearzero(err[i-2], err[i-1], err[i])
      # at a turning point, so potentially two roots between i-2 and i, try and find both...
      lims = _ordered(θ[i-2], θ[i])
      soln = solve(NonlinearProblem{false}(_Δz, T1(θ[i-1]), (pm, tx, p2.x, p2.z)); abstol=pm.atol)
      if successful_retcode(soln.retcode) && abs(soln.resid) < ztol && lims[1] < soln.u < lims[2]
        θ₁ = soln.u
        put!(erays, _trace(pm, tx, θ₁, p2.x, pm.ds; paths)[1])
        if θ[i-1] < θ₁
          soln = solve(NonlinearProblem{false}(_Δz, T1(lims[2]), (pm, tx, p2.x, p2.z)); abstol=pm.atol)
          if successful_retcode(soln.retcode) && abs(soln.resid) < ztol && θ₁ < soln.u < lims[2]
            put!(erays, _trace(pm, tx, soln.u, p2.x, pm.ds; paths)[1])
          end
        else
          soln = solve(NonlinearProblem{false}(_Δz, T1(lims[1]), (pm, tx, p2.x, p2.z)); abstol=pm.atol)
          if successful_retcode(soln.retcode) && abs(soln.resid) < ztol && lims[1] < soln.u < θ₁
            put!(erays, _trace(pm, tx, soln.u, p2.x, pm.ds; paths)[1])
          end
        end
      end
    end
  end
  close(erays)
  sort(collect(erays); by=Base.Fix2(getfield, :t))
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
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  δθ = Float64(θ.step)
  f = frequency(tx)
  h = maximum(pm.env.bathymetry)
  mid_temp = (minimum(pm.env.temperature) + maximum(pm.env.temperature)) / 2
  c₀ = value(pm.env.soundspeed, location(tx))
  rmax = maximum(rxs.xrange) + 0.1
  γ = absorption(f, 1, pm.env.salinity, mid_temp, h/2)  # nominal absorption per meter
  log_γ = log(γ)
  T = promote_type(env_type(pm.env), eltype(location(tx)), typeof(f), ComplexF64)
  afld = zeros(T, size(rxs,1), size(rxs,2))
  afld_lock = [ReentrantLock() for _ ∈ 1:size(rxs,1)]
  Threads.@threads for θ₀ ∈ θ
    ray, aux_info = _trace(pm, tx, θ₀, rmax, pm.ds; aux=true)
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
      C2 = δθ * cos(θ₀) * cₛ / c₀
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
        n = norm(cpa - rxpos)                           # normal distance from ray to rx
        n < 4W || continue                              # rx too far from ray segment
        A = C1 * sqrt(C2 / (cpa[1] * W))                # [COA (3.76)]
        if mode === :coherent
          t = t1 + α * (t2 - t1)                        # time at cpa
          P = A * exp(-(n / W)^2) * cispi(2f * t)       # [COA (3.72), COA (3.75)]
        else
          P = complex(abs2(A * exp(-(n / W)^2)), 0.0)
        end
        lock(afld_lock[j])
        try
          afld[j,k] += P
        finally
          unlock(afld_lock[j])
        end
      end
    end
  end
  mode === :incoherent && (afld .= sqrt.(afld))
  afld * db2amp(spl(tx))
end

### helper functions

_ordered(a, b) = a < b ? (a, b) : (b, a)

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

function _∂u((r, z, ξ, ζ, t, p, q), (c, ∂c, ∂²c), s)
  # implementation based on [COA (3.161-164, 3.58-63)]
  # assumes range-independent soundspeed
  cᵥ = c(z)
  cᵥ² = cᵥ * cᵥ
  c̄ = ∂²c(z) * ξ * ξ
  SA[cᵥ * ξ, cᵥ * ζ, 0, -∂c(z) / cᵥ², 1 / cᵥ, -c̄ * q, cᵥ * p]
end

function _check_ray!(out, u, s, integrator, a, b, rmax)
  out[1] = value(a, (u[1], 0.0, 0.0)) - u[2]    # surface reflection
  out[2] = u[2] + value(b, (u[1], 0.0, 0.0))    # bottom reflection
  out[3] = rmax - u[1]                          # maximum range
  out[4] = u[3]                                 # ray turned back
end

# trace a ray starting at (r0, z0) with angle θ until it hits the surface,
# bottom, or reaches rmax. T is the type to use for computations, and p0, q0
# are the initial spreading parameters (default to 1/c₀ and 0, respectively).
# ds controls the spacing of points along the ray (0 = only events).
# Returns a 2-tuple of vectors containing distances along the ray, and
# (r, z, ξ, ζ, t, p, q) values at each point.
function _trace_segment(T, pm, r0, z0, θ, rmax, ds, p0, q0)
  a = pm.env.altimetry
  b = pm.env.bathymetry
  c = z -> value(pm.env.soundspeed, z)
  ∂c = z -> derivative(c, z)
  ∂²c = z -> derivative(∂c, z)
  cᵥ = c(z0)
  u0 = SA[convert(T, r0), z0, cos(θ)/cᵥ, sin(θ)/cᵥ, zero(T), p0, q0]
  prob = ODEProblem{false}(_∂u, u0, (zero(T), one(T) * pm.rugosity * (rmax-r0)/cos(θ)), (c, ∂c, ∂²c))
  cb = VectorContinuousCallback(
    (out, u, s, i) -> _check_ray!(out, u, s, i, a, b, rmax),
    (i, ndx) -> terminate!(i), 4;
    rootfind = true
  )
  if ds ≤ 0
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, save_everystep=false, callback=cb)
  else
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, saveat=ds, callback=cb)
  end
  (soln.t, soln.u)
end

# trace a ray starting at tx1 with angle θ until it reaches rmax. ds controls the
# spacing of points along the ray. If paths=true, also returns the ray path as a
# vector of (x, y, z) points. If aux=true, also returns auxiliary info as a vector
# of (s, q, t, A) tuples at each point.
function _trace(pm::RaySolver, tx1::AbstractAcousticSource, θ, rmax, ds=0.0; paths=true, aux=false)
  θ₀ = θ
  f = frequency(tx1)
  mid_temp = (minimum(pm.env.temperature) + maximum(pm.env.temperature)) / 2
  zmin = -maximum(pm.env.bathymetry)
  ϵ = one(zmin) * pm.solver_tol
  p = location(tx1)
  c₀ = value(pm.env.soundspeed, p)
  T = promote_type(env_type(pm.env), eltype(p), typeof(f), typeof(θ), typeof(rmax))
  raypath = Array{@NamedTuple{x::T,y::T,z::T}}(undef, 0)
  aux_info = Array{Tuple{T,T,T,Complex{T}}}(undef, 0)
  A = one(Complex{T})   # phasor
  s = 0                 # surface bounces
  b = 0                 # bottom bounces
  t = zero(T)           # time along ray
  D = zero(T)           # distance along ray
  q = zero(T)           # spreading factor
  qp = one(T) / c₀      # spreading rate
  while true
    svec, u = _trace_segment(T, pm, p[1], p[3], θ, rmax, ds, qp, q)
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
    t += dt
    D += dD
    A *= cis(-π/2 * kmah)                 # [COA §3.4.1 (KMAH correction)]
    zmin = -value(pm.env.bathymetry, (r, 0.0, 0.0))
    zmax = value(pm.env.altimetry, (r, 0.0, 0.0))
    p = (r, 0.0, clamp(z, zmin+ϵ, zmax-ϵ))
    if r ≥ rmax - 1e-3
      break
    elseif isapprox(z, zmax; atol=1e-3)   # hit the surface
      s += 1
      ρ = value(pm.env.density, (r, 0.0, 0.0))
      c = value(pm.env.soundspeed, (r, 0.0, 0.0))
      A *= reflection_coef(pm.env.surface, f, π/2 - θ, ρ, c)
      α = atan(derivative(x -> value(pm.env.altimetry, (x, 0.0, 0.0)), r))
      θ = -θ + 2α
    elseif isapprox(z, zmin; atol=1e-3)   # hit the bottom
      b += 1
      ρ = value(pm.env.density, (r, 0.0, z))
      c = value(pm.env.soundspeed, (r, 0.0, z))
      A *= reflection_coef(pm.env.seabed, f, π/2 + θ, ρ, c)
      α = atan(derivative(x -> value(pm.env.bathymetry, (x, 0.0, 0.0)), r))
      θ = -θ - 2α
    else
      break
    end
    abs(θ) ≥ π/2 && break
    abs(A)/abs(q) < pm.min_amplitude && break
  end
  cₛ = value(pm.env.soundspeed, p)
  A *= √abs(cₛ * cos(θ₀) / (p[1] * c₀ * q))                   # [COA (3.65)]
  A *= absorption(f, D, pm.env.salinity, mid_temp, -zmin/2)   # nominal absorption
  RayArrival(t, A, s, b, θ₀, -θ, raypath), aux_info
end

_Δz(ϕ, (pm, tx1, rmax, z)) = _trace(pm, tx1, ϕ, rmax)[1].path[end].z - z
