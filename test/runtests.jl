using AcousticRayTracers
using AcousticRayTracers.UnderwaterAcoustics
using Test
using ForwardDiff

@testset "pm-raysolver" begin

    @test RaySolver in models()

    env = UnderwaterEnvironment()
    pm = RaySolver(env)
    @test pm isa RaySolver

    arr = arrivals(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    @test arr isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    @test length(arr) >= 7
    @test arr[1].time ≈ 0.0650 atol=0.0001
    @test arr[2].time ≈ 0.0657 atol=0.0001
    @test arr[3].time ≈ 0.0670 atol=0.0001
    @test all([arr[j].time > arr[j-1].time for j ∈ 2:7])
    @test abs(arr[1].phasor) ≈ 0.01 atol=0.001
    @test real(arr[2].phasor) < 0.0
    @test imag(arr[2].phasor) ≈ 0.0
    @test all([abs(arr[j].phasor) < abs(arr[j-1].phasor) for j ∈ 2:7])
    @test [(arr[j].surface, arr[j].bottom) for j ∈ 1:7] == [(0,0), (1,0), (0,1), (1,1), (1,1), (2,1), (1,2)]
    @test abs.([a.arrivalangle for a ∈ arr]) ≈ abs.([a.launchangle for a ∈ arr])

    r = eigenrays(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    @test length(r) >= 7
    @test abs.([a.arrivalangle for a ∈ r]) ≈ abs.([a.launchangle for a ∈ r])
    @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:7])
    #@test all([r[j].raypath[end][k] .≈ (100.0, 0.0, -10.0)[k] for j ∈ 1:7, k ∈ 1:3])

    r = rays(pm, AcousticSource(0.0, -5.0, 1000.0), -60°:15°:60°, 100.0)
    @test r isa AbstractArray{<:UnderwaterAcoustics.RayArrival}
    @test length(r) == 9
    @test all([r[j].launchangle for j ∈ 1:9] .≈ -60°:15°:60°)
    @test abs.([a.arrivalangle for a ∈ r]) ≈ abs.([a.launchangle for a ∈ r])
    @test all([r[j].raypath[1] == (0.0, 0.0, -5.0) for j ∈ 1:9])
    @test r[4].raypath[end][1] ≥ 99.9
    @test r[5].raypath[end][1] ≥ 99.9
    @test r[6].raypath[end][1] ≥ 99.9
    @test r[7].raypath[end][1] ≥ 99.9

    env = UnderwaterEnvironment(ssp=IsoSSP(1500.0), seabed=RayleighReflectionCoef(1.0, 1.0))
    pm = RaySolver(env)
    d = (√1209.0)/4.0
    x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
    @test x isa Complex
    @test abs(x) ≈ 0.0 atol=0.0002
    x′ = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
    @test x′ isa Complex
    @test imag(x′) == 0.0
    @test abs(x′) > 1/100.0
    d = (√2409.0)/8.0
    x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
    @test abs(x) > abs(x′)
    y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d))
    @test -10 * log10(abs2(x)) ≈ y atol=0.1
    x = transfercoef(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
    @test abs(x) ≈ abs(x′) atol=0.0001
    y = transmissionloss(pm, AcousticSource(0.0, -d, 1000.0), AcousticReceiver(100.0, -d); mode=:incoherent)
    @test -10 * log10(abs2(x)) ≈ y atol=0.1
    x1 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    x2 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    x3 = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1 x2 x3] ≈ x atol=0.01
    y = transfercoef(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
    @test y isa AbstractMatrix
    @test size(y) == (3, 3)
    @test x' ≈ y[1,:] atol=0.05
    x1 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -5.0))
    x2 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -10.0))
    x3 = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiver(100.0, -15.0))
    x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), [AcousticReceiver(100.0, -d) for d ∈ 5.0:5.0:15.0])
    @test x isa AbstractVector
    @test [x1, x2, x3] == x
    x = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 0.0, 1, -5.0, -5.0, 3))
    @test x isa AbstractMatrix
    @test size(x) == (1, 3)
    @test [x1, x2] ≈ x[1:2] atol=1.5
    y = transmissionloss(pm, AcousticSource(0.0, -5.0, 1000.0), AcousticReceiverGrid2D(100.0, 10.0, 3, -5.0, -5.0, 3))
    @test y isa AbstractMatrix
    @test size(y) == (3, 3)
    @test x' ≈ y[1,:] atol=0.1

    tx = AcousticSource(0.0, -5.0, 1000.0)
    rx = AcousticReceiver(100.0, -10.0)
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
    rx = AcousticReceiver(100.0, -10.0)
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,)
    tx = [AcousticSource(0.0, -5.0, 1000.0), AcousticSource(0.0, -10.0, 2000.0)]
    rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    tx = AcousticSource(0.0, -5.0, 1000.0)
    rx = [AcousticReceiver(100.0, -10.0), AcousticReceiver(100.0, -15.0)]
    sig = record(pm, tx, rx, 1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = record(pm, tx, rx, 1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0)
    @test size(sig) == (44100,2)
    sig = recorder(pm, tx, rx)(1.0, 44100.0; start=0.5)
    @test size(sig) == (44100,2)

  end

  function ∂(f, x, i, ϵ)
    x1 = copy(x)
    x1[i] = x[i] + ϵ
    f1 = f(x1)
    x1[i] = x[i] - ϵ
    (f1 - f(x1)) / (2ϵ)
  end

  @testset "pm-∂raysolver" begin

    function ℳ₁(x)
      D, R, d1, d2, f, c = x
      env = UnderwaterEnvironment(ssp=IsoSSP(c), bathymetry=ConstantDepth(D))
      pm = RaySolver(env)
      transmissionloss(pm, AcousticSource(0.0, -d1, f), AcousticReceiver(R, -d2))
    end

    function ℳ₂(x)
      D, R, d1, d2, f, c = x
      env = UnderwaterEnvironment(ssp=IsoSSP(c), bathymetry=ConstantDepth(D))
      pm = RaySolver(env)
      transmissionloss(pm, AcousticSource(0.0, -d1, f), AcousticReceiverGrid2D(R, 0.0, 1, -d2, 0.0, 1))[1,1]
    end

    x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
    ∇ℳ = ForwardDiff.gradient(ℳ₁, x)
    for i ∈ 1:length(x)
      # skip i = 2 because it is not yet supported
      i != 2 && @test ∇ℳ[i] ≈ ∂(ℳ₁, x, i, 0.0001) atol=0.1
    end

    x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
    ∇ℳ = ForwardDiff.gradient(ℳ₁, x)
    for i ∈ 1:length(x)
      # skip i = 2 because it is not yet supported
      i != 2 && @test ∇ℳ[i] ≈ ∂(ℳ₁, x, i, 0.0001) atol=0.1
    end

    x = [20.0, 100.0, 5.0, 10.0, 5000.0, 1500.0]
    ∇ℳ = ForwardDiff.gradient(ℳ₂, x)
    for i ∈ 1:length(x)
      @test ∇ℳ[i] ≈ ∂(ℳ₂, x, i, 0.0001) atol=0.1
    end

    x = [25.0, 200.0, 10.0, 8.0, 1000.0, 1540.0]
    ∇ℳ = ForwardDiff.gradient(ℳ₂, x)
    for i ∈ 1:length(x)
      @test ∇ℳ[i] ≈ ∂(ℳ₂, x, i, 0.0001) atol=0.1
    end

  end
