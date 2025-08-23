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
  ssp = SampledField([
    1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1476.7, 1472.6, 1468.8,
    1467.2, 1471.6, 1473.6, 1473.6, 1472.7, 1472.2, 1471.6, 1471.6, 1472.0, 1472.7, 1473.1,
    1474.9, 1477.0, 1478.1, 1480.7, 1483.8, 1490.5, 1498.3, 1506.5]; z=[0, -5, -10, -15, -20,
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
  @test xloss[50,50] ≈ 71.8 atol=0.5
  @test xloss[100,100] ≈ 73.4 atol=0.5
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
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂
  x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
  ∇ℳ₁ = gradient(ℳ₁, fd, x)
  ∇ℳ₂ = gradient(ℳ₂, fd, x)
  @test gradient(ℳ₁, AutoForwardDiff(), x) ≈ ∇ℳ₁
  @test gradient(ℳ₂, AutoForwardDiff(), x) ≈ ∇ℳ₂
end
