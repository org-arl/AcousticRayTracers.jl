using TestItems

@testitem "models" begin
  using UnderwaterAcoustics
  m = models()
  @test m isa Vector{Type{<:UnderwaterAcoustics.AbstractPropagationModel}}
  @test RaySolver ∈ m
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", temperature = 27.0u"°C", salinity = 35.0u"ppt", seabed = SandySilt)
end

@testitem "raysolver-arrivals+ir" begin
  using UnderwaterAcoustics
  using UnderwaterAcoustics: distance
  env = UnderwaterEnvironment(bathymetry=20)
  pm = @inferred RaySolver(env)
  @test pm isa RaySolver
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = AcousticReceiver(100.0, -10.0)
  arr = @inferred arrivals(pm, tx, rx)
  @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
  @test length(arr) ≥ 7
  @test arr[1].t ≈ 0.0650 atol=0.0001
  @test arr[2].t ≈ 0.0657 atol=0.0001
  @test arr[3].t ≈ 0.0670 atol=0.0001
  @test all([arr[j].t > arr[j-1].t for j ∈ 2:7])
  @test abs(arr[1].ϕ) ≈ 0.01 atol=0.001
  @test real(arr[2].ϕ) < 0.0
  @test imag(arr[2].ϕ) ≈ 0.0
  @test all([abs(arr[j].ϕ) < abs(arr[j-1].ϕ) for j ∈ 2:7])
  @test [(arr[j].ns, arr[j].nb) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
  @test abs.([a.θᵣ for a ∈ arr]) ≈ abs.([a.θₛ for a ∈ arr])
  @test all([a.path[1] == (x=0.0, y=0.0, z=-5.0) for a ∈ arr])
  @test all([distance(a.path[end], (x=100.0, y=0.0, z=-10.0)) < 0.5 for a ∈ arr])
  ir1 = @inferred impulse_response(pm, tx, rx, 10000.0; abstime=false)
  ir2 = @inferred impulse_response(pm, tx, rx, 10000.0; abstime=true)
  @test ir1 isa AbstractVector{<:Complex}
  @test ir2 isa AbstractVector{<:Complex}
end

@testitem "raysolver-field" begin
  using UnderwaterAcoustics
  env = UnderwaterEnvironment(bathymetry = 20.0u"m", soundspeed = 1500.0u"m/s",
    seabed = FluidBoundary(water_density(), 1500.0u"m/s", 0.0))  # non-reflecting
  pm = @inferred RaySolver(env)
  d = (√1209.0)/4.0
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test x isa Complex
  @test abs(x) ≈ 0.0 atol=0.0002
  x′ = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test x′ isa Complex
  @test imag(x′) == 0.0
  @test abs(x′) > 1/100.0
  d = (√2409.0)/8.0
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test abs(x) > abs(x′)
  y = @inferred transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)))
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  x = @inferred acoustic_field(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test abs(x) ≈ abs(x′) atol=0.0001
  y = @inferred transmission_loss(pm, AcousticSource((x=0.0, z=-d), 1000.0), AcousticReceiver((x=100.0, z=-d)); mode=:incoherent)
  @test -10 * log10(abs2(x)) ≈ y atol=0.1
  tx = @inferred AcousticSource((x=0.0, z=-5.0), 1000.0)
  x1 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-5.0)))
  x2 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-10.0)))
  x3 = @inferred acoustic_field(pm, tx, AcousticReceiver((x=100.0, z=-15.0)))
  x = @inferred acoustic_field(pm, tx, [AcousticReceiver((x=100.0, z=-d)) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = @inferred acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test [x1 x2 x3] ≈ x atol=0.1
  x = @inferred acoustic_field(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test [x1, x2, x3] ≈ x[1,:] atol=0.1
  x1 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -5.0))
  x2 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -10.0))
  x3 = @inferred transmission_loss(pm, tx, AcousticReceiver(100.0, -15.0))
  x = @inferred transmission_loss(pm, tx, [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
  @test x isa AbstractVector
  @test [x1, x2, x3] == x
  x = @inferred transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:100.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (1, 3)
  @test all(abs.([x1 x2 x3] - x) .< 1.5)
  x = @inferred transmission_loss(pm, tx, AcousticReceiverGrid2D(100.0:10.0:120.0, -5.0:-5.0:-15.0))
  @test x isa AbstractMatrix
  @test size(x) == (3, 3)
  @test all(abs.([x1, x2, x3] - x[1,:]) .< 1.5)
end

@testitem "raysolver-seamount-tl" begin
  using UnderwaterAcoustics
  bathy = SampledField([3000, 3000, 500, 3000, 3000]; x=[0, 10e3, 20e3, 30e3, 100e3])
  # FIXME although ssp should start at z=0, it starts at z=1e-6 because evaluation
  # of derivative in RaySolver at z=0 in the former case seems to cause an
  # exception only inside a @testitem, but causes no problems in a REPL or a script
  # Upstream bug: https://github.com/JuliaMath/Interpolations.jl/issues/645
  ssp = SampledField([
    1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1472.6, 1468.8,
    1467.2, 1471.6, 1473.6, 1473.6, 1472.7, 1472.2, 1471.6, 1471.6, 1472.0, 1472.7, 1473.1,
    1474.9, 1477.0, 1478.1, 1480.7, 1483.8, 1490.5, 1498.3, 1506.5]; z=[1e-6, -5, -10, -15, -20,
    -25, -30, -35, -38, -50, -70, -100, -140, -160, -170, -200, -215, -250, -300, -370, -450,
    -500, -700, -900, -1000, -1250, -1500, -2000, -2500, -3000])
  env = UnderwaterEnvironment(
    bathymetry = bathy,
    soundspeed = ssp,
    seabed = FluidBoundary(1.5*water_density(), 1550.0, dBperλ(0.5)))
  tx = AcousticSource((x=0, z=-18), 230.0)
  rxs = AcousticReceiverGrid2D(0:100:100000, -3000:15:0)
  pm = @inferred RaySolver(env; min_angle=-89°, max_angle=89°)
  xloss = @inferred transmission_loss(pm, tx, rxs; mode=:incoherent)
  @test size(xloss) == (1001, 201)
  @test xloss[200,50] > 150
  # reference values at the auto-nbeams default (cap raised 1000 → 4000; the
  # old constants 71.1 / 76.6 were unconverged-in-beams artifacts); BELLHOP
  # gives 75.69 at [100,100] — the remaining bias is the SSP-knot p-jump
  # (issue #24) and closes to 75.77 with its fix
  @test xloss[50,50] ≈ 69.98 atol=0.5
  @test xloss[100,100] ≈ 73.48 atol=0.5
end

@testitem "raysolver-backscatter" begin
  using UnderwaterAcoustics
  # null test: with nothing to backscatter, arrivals match the default solver
  env = UnderwaterEnvironment(bathymetry=20)
  tx = AcousticSource(0.0, -5.0, 1000.0)
  rx = AcousticReceiver(100.0, -10.0)
  arr0 = arrivals(RaySolver(env), tx, rx)
  pm = @inferred RaySolver(env; backscatter=true)
  arr1 = arrivals(pm, tx, rx)
  @test length(arr1) ≥ 7
  for a0 ∈ arr0[1:7]
    @test any(a -> isapprox(a.t, a0.t; atol=1e-4) && isapprox(abs(a.ϕ), abs(a0.ϕ); rtol=0.05), arr1)
  end
  @test all(a -> abs(a.θᵣ) < π/2, arr1)
  # bathymetric "wall" behind the receiver: rays travel past the receiver,
  # reverse off the rigid slope (slope + surface bounce combinations) and
  # return as backscattered arrivals (|arrival angle| > π/2)
  bathy = SampledField([20, 20, 1]; x=[0, 85, 100])
  env2 = UnderwaterEnvironment(bathymetry=bathy, seabed=RigidBoundary, soundspeed=1500.0)
  tx2 = AcousticSource(0.0, -10.0, 1000.0)
  rx2 = AcousticReceiver(40.0, -10.0)
  pm2 = RaySolver(env2; backscatter=true, rmax=100.0)
  arr2 = arrivals(pm2, tx2, rx2)
  fwd = sort(filter(a -> abs(a.θᵣ) ≤ π/2, arr2); by=a->a.t)
  bck = sort(filter(a -> abs(a.θᵣ) > π/2, arr2); by=a->a.t)
  @test length(fwd) ≥ 3
  @test length(bck) ≥ 4
  @test fwd[1].t ≈ 40/1500 atol=1e-4                # direct arrival
  @test bck[1].t ≈ 0.1106 atol=1e-3                 # earliest wall echo
  @test bck[1].t > (2*85 - 40)/1500                 # slower than a straight retro-reflection
  @test abs(bck[1].ϕ) ≈ 0.0114 rtol=0.1             # echo amplitude
  @test abs(bck[1].ϕ) < abs(fwd[1].ϕ)               # weaker than the direct arrival
  # all eigenray paths (forward and backscattered) must end at the receiver
  @test all(a -> abs(a.path[end].x - 40) < 0.5 && abs(a.path[end].z + 10) < 0.5, arr2)
  # with default rmax (= query range), the wall at 85-100 m is outside the
  # domain and produces no echoes
  arr3 = arrivals(RaySolver(env2; backscatter=true), tx2, rx2)
  @test count(a -> abs(a.θᵣ) > π/2, arr3) == 0
  # wall to the LEFT of the source: echoes require rmin < 0 (mirrored fan)
  bathy4 = SampledField([2, 20, 20]; x=[-100, -80, 200])
  env4 = UnderwaterEnvironment(bathymetry=bathy4, seabed=RigidBoundary, soundspeed=1500.0)
  tx4 = AcousticSource(0.0, -10.0, 1000.0)
  rx4 = AcousticReceiver(50.0, -10.0)
  arr4a = arrivals(RaySolver(env4; backscatter=true), tx4, rx4)
  arr4b = arrivals(RaySolver(env4; backscatter=true, rmin=-100.0, rmax=100.0), tx4, rx4)
  @test count(a -> abs(a.θₛ) > π/2, arr4a) == 0
  @test count(a -> abs(a.θₛ) > π/2, arr4b) > 0
  @test all(a -> abs(a.path[end].x - 50) < 0.5 && abs(a.path[end].z + 10) < 0.5, arr4b)
end

@testitem "∂raysolver-backscatter" begin
  using UnderwaterAcoustics
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = RaySolver(env; backscatter=true)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
  @test gradient(ℳ, AutoForwardDiff(), x) ≈ gradient(ℳ, fd, x) atol=1e-3
end

@testitem "raysolver-scatterers" begin
  using UnderwaterAcoustics
  if !isdefined(UnderwaterAcoustics, :Scatterer)
    # awaiting the UnderwaterAcoustics.jl scatterer API (rev 3)
    @test_skip false
  else
    # backscatter off a rigid circular scatterer: expect an arrival at
    # t ≈ (48 + 38) / 1500 (tx at 0, rx at 10, scatterer surface at 48)
    env = UnderwaterEnvironment(bathymetry = 40.0, soundspeed = 1500.0,
      scatterers = (Scatterer(Ellipse(50.0, -20.0, 2.0, 2.0), RigidBoundary),))
    tx = AcousticSource(0.0, -20.0, 5000.0)
    rx = AcousticReceiver(10.0, -20.0)
    pm = RaySolver(env; backscatter=true, rmax=60.0)
    arr = arrivals(pm, tx, rx)
    bsc = filter(a -> abs(a.θᵣ) > π/2, arr)
    @test any(a -> isapprox(a.t, 86/1500; atol=1e-3), bsc)
    # pressure-release scatterer: backscattered arrival flips phase by π
    env2 = UnderwaterEnvironment(bathymetry = 40.0, soundspeed = 1500.0,
      scatterers = (Scatterer(Ellipse(50.0, -20.0, 2.0, 2.0), PressureReleaseBoundary),))
    arr2 = arrivals(RaySolver(env2; backscatter=true, rmax=60.0), tx, rx)
    bsc2 = filter(a -> abs(a.θᵣ) > π/2, arr2)
    a1 = bsc[argmin([abs(a.t - 86/1500) for a ∈ bsc])]
    a2 = bsc2[argmin([abs(a.t - 86/1500) for a ∈ bsc2])]
    @test abs(angle(a1.ϕ / a2.ϕ)) ≈ π atol=0.1
    # anti-tunneling: a small scatterer (0.5 m radius) must still be detected
    env3 = UnderwaterEnvironment(bathymetry = 40.0, soundspeed = 1500.0,
      scatterers = (Scatterer(Ellipse(50.0, -20.0, 0.5, 0.5), RigidBoundary),))
    arr3 = arrivals(RaySolver(env3; backscatter=true, rmax=60.0), tx, rx)
    @test any(a -> abs(a.θᵣ) > π/2 && isapprox(a.t, (49.5 + 39.5)/1500; atol=1e-3), arr3)
    # shadow: coherent field behind the scatterer drops relative to no-scatterer case
    env4 = UnderwaterEnvironment(bathymetry = 40.0, soundspeed = 1500.0)
    rx2 = AcousticReceiver(60.0, -20.0)
    tl_free = transmission_loss(RaySolver(env4), tx, rx2)
    tl_shad = transmission_loss(RaySolver(env; backscatter=true, rmax=60.0), tx, rx2)
    @test tl_shad > tl_free
    # AD through scatterer geometry
    using DifferentiationInterface
    import ForwardDiff, FiniteDifferences
    fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
    function ℳ((sx, sz, sa))
      envg = UnderwaterEnvironment(bathymetry = 40.0, soundspeed = 1500.0,
        scatterers = (Scatterer(Ellipse(sx, sz, sa, sa), RigidBoundary),))
      pmg = RaySolver(envg; backscatter=true, rmax=60.0)
      transmission_loss(pmg, AcousticSource((x=0.0, z=-20.0), 5000.0), AcousticReceiver((x=10.0, z=-20.0)))
    end
    x = [50.0, -20.0, 2.0]
    @test gradient(ℳ, AutoForwardDiff(), x) ≈ gradient(ℳ, fd, x) atol=1e-3
  end
end

@testitem "raysolver-curvature-benchmark" begin
  using UnderwaterAcoustics
  using UnderwaterAcoustics: absorption
  if !isdefined(UnderwaterAcoustics, :Scatterer)
    # awaiting the UnderwaterAcoustics.jl scatterer API (rev 3)
    @test_skip false
  else
    # verify the curved-boundary beam-spreading correction against the mirror
    # equation (geometric optics, independent of the solver's dynamic ray
    # tracing): a convex circular reflector transforms the in-plane wavefront
    # curvature as 1/ρ' = 1/s₁ + 2κ/sinθg
    c = 1500.0; f = 5000.0; a = 2.0; L = 50.0
    env = UnderwaterEnvironment(bathymetry=200.0, soundspeed=c,
      seabed=FluidBoundary(water_density(), c, 0.0),  # non-reflecting
      scatterers=(Scatterer(Ellipse(L, -100.0, a, a), RigidBoundary),))
    mid_temp = (minimum(env.temperature) + maximum(env.temperature)) / 2
    hmax = 200.0
    tx = AcousticSource(0.0, -100.0, f)
    pm = RaySolver(env; backscatter=true, rmax=100.0, nbeams=2000)
    function analytic_echo(txp, rxp)
      # specular point on the tx-facing half of the circle by Fermat
      path(u) = begin
        p = boundary_point(Ellipse(L, -100.0, a, a), u)
        hypot(p.x - txp[1], p.z - txp[2]) + hypot(p.x - rxp[1], p.z - rxp[2])
      end
      lo, hi = 0.25, 0.75
      φ = (√5 - 1) / 2
      for _ in 1:200
        m1 = hi - φ * (hi - lo); m2 = lo + φ * (hi - lo)
        path(m1) < path(m2) ? (hi = m2) : (lo = m1)
      end
      u = (lo + hi) / 2
      P = boundary_point(Ellipse(L, -100.0, a, a), u)
      n̂ = (x=(P.x - L)/a, z=(P.z + 100.0)/a)
      s₁ = hypot(P.x - txp[1], P.z - txp[2])
      s₂ = hypot(P.x - rxp[1], P.z - rxp[2])
      t̂ = ((P.x - txp[1])/s₁, (P.z - txp[2])/s₁)
      sinθg = abs(t̂[1]*n̂.x + t̂[2]*n̂.z)
      ρ′ = 1 / (1/s₁ + 2/(a*sinθg))          # mirror equation, convex
      q_rx = s₁ * (ρ′ + s₂) / ρ′             # in-plane spreading at receiver
      θ₀ = atan(P.z - txp[2], P.x - txp[1])
      D = s₁ + s₂
      A = √abs(cos(θ₀) * c / (rxp[1] * c * q_rx)) * absorption(f, D, env.salinity, mid_temp, hmax/2)
      (A=A, t=D/c)
    end
    # normal incidence through strongly oblique (sinθg ≈ 0.93)
    for rxp ∈ ((10.0, -100.0), (15.0, -92.0), (10.0, -70.0), (5.0, -60.0))
      arr = arrivals(pm, tx, AcousticReceiver(rxp[1], rxp[2]))
      ref = analytic_echo((0.0, -100.0), rxp)
      bck = filter(x -> abs(x.θᵣ) > π/2 && isapprox(x.t, ref.t; atol=2e-4), arr)
      @test !isempty(bck)
      if !isempty(bck)
        m = bck[argmin([abs(x.t - ref.t) for x ∈ bck])]
        @test abs(m.ϕ) ≈ ref.A rtol=1e-3
      end
    end
  end
end

@testitem "∂raysolver" begin
  using UnderwaterAcoustics
  using DifferentiationInterface
  import ForwardDiff, FiniteDifferences
  fd = AutoFiniteDifferences(fdm=FiniteDifferences.central_fdm(5, 1))
  function ℳ₁((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = RaySolver(env)
    transmission_loss(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)))
  end
  function ℳ₂((D, R, d1, d2, f, c))
    env = UnderwaterEnvironment(bathymetry = D, soundspeed = c, seabed = SandySilt)
    pm = RaySolver(env)
    ir = impulse_response(pm, AcousticSource((x=0.0, z=-d1), f), AcousticReceiver((x=R, z=-d2)), 10000.0; abstime=true)
    sum(abs2, samples(ir))
  end
  x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
  ∇ℳ₁ = gradient(ℳ₁, fd, x)
  ∇ℳ₂ = gradient(ℳ₂, fd, x)
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁ atol=1e-3
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂ atol=1e-3
  x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
  ∇ℳ₁ = gradient(ℳ₁, fd, x)
  ∇ℳ₂ = gradient(ℳ₂, fd, x)
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁ atol=1e-3
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂ atol=1e-3
end
