import LinearAlgebra: norm, dot
import OrdinaryDiffEq: ODEProblem, VectorContinuousCallback, Tsit5, Rodas5, solve, terminate!
import ForwardDiff: derivative
import NonlinearSolve: IntervalNonlinearProblem
import NonlinearSolve.SciMLBase: successful_retcode

Base.@kwdef struct RaySolver{T1,T2} <: AbstractRayPropagationModel
  env::T1
  nbeams::Int = 0
  min_angle::Float64 = -deg2rad(80)
  max_angle::Float64 = +deg2rad(80)
  ds::Float64 = 1.0
  atol::Float64 = 1e-4
  rugosity::Float64 = 1.5
  min_amplitude::Float64 = 1e-5
  solver::T2 = nothing
  solver_tol::Float64 = 1e-4
  function RaySolver(env, nbeams, min_angle, max_angle, ds, atol, rugosity, min_amplitude, solver, solver_tol)
    nbeams < 0 && (nbeams = 0)
    -π/2 ≤ min_angle ≤ π/2 || error("min_angle should be between -π/2 and π/2")
    -π/2 ≤ max_angle ≤ π/2 || error("max_angle should be between -π/2 and π/2")
    min_angle < max_angle || error("max_angle should be more than min_angle")
    solver = something(solver, is_isovelocity(env) ? Tsit5() : Rodas5())
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
- `rugosity`: rugosity of the seabed (default: 1.5 m)
- `min_amplitude`: minimum ray amplitude to track (default: 1e-5)
- `solver`: differential equation solver (default: nothing, auto)
- `solver_tol`: solver tolerance (default: 1e-4)
"""
RaySolver(env; kwargs...) = RaySolver(; env=env, kwargs...)

Base.show(io::IO, pm::RaySolver) = print(io, "RaySolver(⋯)")

### interface functions

function UnderwaterAcoustics.arrivals(pm::RaySolver, tx::AbstractAcousticSource, rx::AbstractAcousticReceiver; paths=true)
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
  rays = tmap(θ1 -> _trace(pm, tx, θ1, p2.x), θ)
  err = map(rays) do ray
    p3 = ray.path[end]
    isapprox(p3.x, p2.x; atol=pm.atol) && isapprox(p3.y, p2.y; atol=pm.atol) ? p3.z - p2.z : NaN
  end
  T1 = promote_type(env_type(pm.env), eltype(location(tx)), typeof(p2.x), typeof(p2.z))
  T2 = T1 == eltype(θ) ? eltype(rays) : typeof(_trace(pm, tx, T1(θ[1]), p2.x; paths))
  erays = T2[]
  for i ∈ 1:length(θ)
    if isapprox(err[i], 0; atol=pm.atol)
      push!(erays, T1 == eltype(θ) ? rays[i] : _trace(pm, tx, T1(θ[i]), p2.x; paths))
    elseif i > 1 && !isnan(err[i-1]) && !isnan(err[i]) && sign(err[i-1]) * sign(err[i]) < 0
      soln = solve(IntervalNonlinearProblem{false}(_Δz, T1.(_ordered(θ[i-1], θ[i])), (pm, tx, p2.x, p2.z)))
      successful_retcode(soln.retcode) && push!(erays, _trace(pm, tx, soln.u, p2.x, pm.ds; paths))
    # FIXME elseif i > 2 && _isnearzero(err[i-2], err[i-1], err[i])
    #   soln = solve(IntervalNonlinearProblem{false}(_Δz, T1.(_ordered(θ[i-2], θ[i])), (pm, tx, p2.x, p2.z)))
    #   successful_retcode(soln.retcode) && push!(erays, _trace(pm, tx, soln.u, p2.x, pm.ds; paths))
    end
  end
  sort(erays; by=Base.Fix2(getfield, :t))
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

# implementation primarily based on ideas from COA (Computational Ocean Acoustics, 2nd ed., ch. 3)
function UnderwaterAcoustics.acoustic_field(pm::RaySolver, tx::AbstractAcousticSource, rxs::AcousticReceiverGrid2D; mode=:coherent)
  _check2d([tx], rxs)
  mode ∈ (:coherent, :incoherent) || error("Unknown mode :" * string(mode))
  f = frequency(tx)
  min_temp = minimum(pm.env.temperature)
  T = promote_type(env_type(pm.env), eltype(location(tx)), typeof(f), ComplexF64)
  tc = zeros(T, size(rxs,1), size(rxs,2), Threads.nthreads())
  rmax = maximum(rxs.xrange) + 0.1
  h = maximum(pm.env.bathymetry)
  if pm.nbeams > 0
    nbeams = pm.nbeams
  else
    p1 = location(tx)
    R = abs(rmax - p1[1])
    nbeams = clamp(ceil(Int, 20 * (pm.max_angle - pm.min_angle) / atan(h, R)), 100, 1000)
  end
  θ = range(pm.min_angle, pm.max_angle; length=nbeams)
  δθ = Float64(θ.step)
  c₀ = value(pm.env.soundspeed, location(tx))
  ω = 2π * f
  G = 1 / (2π)^(1/4)
  Threads.@threads for θ1 ∈ θ
    β = (mode === :incoherent ? 2 : 1) * cos(θ1) / c₀
    _trace(pm, tx, θ1, rmax, 1.0; cb = (s1, u1, s2, u2, A₀, D₀, t₀, cₛ, kmah1, kmah2) -> begin
      r1, z1, ξ1, ζ1, t1, _, q1 = u1
      r2, z2, ξ2, ζ2, t2, _, q2 = u2
      ndx = findall(r -> r1 ≤ r < r2, rxs.xrange)
      rz1 = (r1, z1)
      vlen = norm((r2-r1, z2-z1))
      rvlen = norm((ξ1, ζ1))
      tᵥ = (ξ1, ζ1) ./ rvlen
      nᵥ = (ζ1, -ξ1) ./ rvlen
      D = D₀ + (s1 + s2) / 2
      γ = absorption(f, D, pm.env.salinity, min_temp, h/2)  # nominal absorption
      Wmax4 = 4 * max(q1, q2) * δθ
      ndx2 = findall(z -> min(z1, z2) - Wmax4 ≤ z ≤ max(z1, z2) + Wmax4, rxs.zrange)
      for j ∈ ndx
        for i ∈ ndx2
          rz = (rxs.xrange[j], rxs.zrange[i])
          v = rz .- rz1
          s = dot(v, tᵥ)
          n = dot(v, nᵥ)
          α = s / vlen
          t = t₀ + t1 + (t2 - t1) * α
          q = q1 + (q2 - q1) * α
          A = A₀ * γ * G * √abs(β * cₛ / (rz[1] * q)) # COA (3.76)
          kmah = round(Int, kmah1 + (kmah2 - kmah1) * α)
          A *= cis(-π/2 * kmah)                       # COA section 3.4.1 (KMAH correction)
          W = abs(q * δθ)                             # COA (3.74)
          A *= exp(-(n / W)^2)
          tc[j, i, Threads.threadid()] += mode === :coherent ? conj(A) * cis(ω * t) : Complex(abs2(A), 0.0)
        end
      end
    end)
  end
  rv = dropdims(sum(tc; dims=3); dims=3)
  mode === :incoherent && (rv = sqrt.(rv))
  rv
end

### helper functions

_ordered(a, b) = a < b ? (a, b) : (b, a)

# function _isnearzero(a, b, c)
#   (isnan(a) || isnan(b) || isnan(c)) && return false
#   sign(a) == sign(b) == sign(c) || return false
#   abs(a) < abs(b) && return false
#   abs(c) < abs(b) && return false
#   return true
# end

function _check2d(tx, rx)
  all(location(tx1).x == 0.0 for tx1 ∈ tx) || error("RaySolver requires transmitters at (0, 0, z)")
  all(location(tx1).y == 0.0 for tx1 ∈ tx) || error("RaySolver requires transmitters in the x-z plane")
  all(location(rx1).x >= 0.0 for rx1 ∈ rx) || error("RaySolver requires receivers to be in the +x halfspace")
  all(location(rx1).y == 0.0 for rx1 ∈ rx) || error("RaySolver requires receivers in the x-z plane")
end

function _∂u!(du, (r, z, ξ, ζ, t, p, q), (c, ∂c, ∂²c), s)
  # implementation based on COA (3.161-164, 3.58-63)
  # assumes range-independent soundspeed
  cᵥ = c(z)
  cᵥ² = cᵥ * cᵥ
  c̄ = ∂²c(z) * ξ * ξ
  du[1] = cᵥ * ξ
  du[2] = cᵥ * ζ
  du[3] = 0
  du[4] = -∂c(z) / cᵥ²
  du[5] = 1 / cᵥ
  du[6] = -c̄ * q
  du[7] = cᵥ * p
end

function _check_ray!(out, u, s, integrator, a, b, rmax)
  out[1] = value(a, (u[1], 0.0, 0.0)) - u[2]    # surface reflection
  out[2] = u[2] + value(b, (u[1], 0.0, 0.0))    # bottom reflection
  out[3] = rmax - u[1]                          # maximum range
  out[4] = u[3]                                 # ray turned back
end

function _trace1(T, pm, r0, z0, θ, rmax, ds, p0, q0)
  a = pm.env.altimetry
  b = pm.env.bathymetry
  c = z -> value(pm.env.soundspeed, z)
  ∂c = z -> derivative(c, z)
  ∂²c = z -> derivative(∂c, z)
  cᵥ = c(z0)
  u0 = [convert(T, r0), z0, cos(θ)/cᵥ, sin(θ)/cᵥ, zero(T), p0, q0]
  prob = ODEProblem{true}(_∂u!, u0, (zero(T), one(T) * pm.rugosity * (rmax-r0)/cos(θ)), (c, ∂c, ∂²c))
  cb = VectorContinuousCallback(
    (out, u, s, i) -> _check_ray!(out, u, s, i, a, b, rmax),
    (i, ndx) -> terminate!(i), 4; rootfind=true)
  if ds ≤ 0
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, save_everystep=false, callback=cb)
  else
    soln = solve(prob, pm.solver; abstol=pm.solver_tol, saveat=ds, callback=cb)
  end
  s2 = soln.u[end]
  soln.t[end], s2[1], s2[2], atan(s2[4], s2[3]), s2[5], s2[6], s2[7], soln.u, soln.t
end

function _trace(pm::RaySolver, tx1::AbstractAcousticSource, θ, rmax, ds=0.0; cb=nothing, paths=true)
  θ₀ = θ
  f = frequency(tx1)
  min_temp = minimum(pm.env.temperature)
  zmin = -maximum(pm.env.bathymetry)
  ϵ = one(zmin) * pm.solver_tol
  p = location(tx1)
  c₀ = value(pm.env.soundspeed, p)
  T = promote_type(env_type(pm.env), eltype(p), typeof(f), typeof(θ), typeof(rmax))
  raypath = Array{typeof(p)}(undef, 0)
  A = one(Complex{T})   # phasor
  s = 0                 # surface bounces
  b = 0                 # bottom bounces
  t = zero(T)           # time along ray
  D = zero(T)           # distance along ray
  q = zero(T)           # spreading factor
  qp = one(T) / c₀      # spreading rate
  while true
    dD, r, z, θ, dt, qp, q, u, svec = _trace1(T, pm, p[1], p[3], θ, rmax, ds, qp, q)
    oq = u[1][7]
    kmah = 0
    kmahhist = zeros(length(u))
    for i ∈ 1:length(u)
      sign(oq) * sign(u[i][7]) == -1 && (kmah += 1)
      kmahhist[i] = kmah
      u[i][7] != 0.0 && (oq = u[i][7])
      paths && push!(raypath, xyz(u[i][1], 0.0, u[i][2]))
    end
    if cb !== nothing
      for i ∈ 2:length(u)
        cₛ = value(pm.env.soundspeed, ((u[i-1][1]+u[i][1])/2, 0.0, (u[i-1][2]+u[i][2])/2))
        cb(svec[i-1], u[i-1], svec[i], u[i], A, D, t, cₛ, kmahhist[i-1], kmahhist[i])
      end
    end
    t += dt
    D += dD
    A *= cis(-π/2 * kmah)                 # COA section 3.4.1 (KMAH correction)
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
    if abs(θ) ≥ π/2
      break
    end
    if abs(A)/abs(q) < pm.min_amplitude
      break
    end
  end
  cₛ = value(pm.env.soundspeed, p)
  A *= √abs(cₛ * cos(θ₀) / (p[1] * c₀ * q))                   # COA (3.65)
  A *= absorption(f, D, pm.env.salinity, min_temp, -zmin/2)   # nominal absorption
  RayArrival(t, A, s, b, θ₀, -θ, raypath)
end

_Δz(ϕ, (pm, tx1, rmax, z)) = _trace(pm, tx1, ϕ, rmax).path[end].z - z
