# Benchmark scenarios: RaySolver vs BELLHOP.
#
# Each scenario is a NamedTuple with:
#   name        — short identifier (used in cache keys and filenames)
#   description — one-line human description
#   tests       — what physics/code-path this scenario probes
#   env         — UnderwaterEnvironment (shared verbatim by both models)
#   tx          — acoustic source
#   grid        — AcousticReceiverGrid2D for TL maps
#   rxs         — a few point receivers for arrivals comparison
#   referee     — :pekeris, :kraken, :lloyd or nothing (independent ground truth)
#
# Scenarios are pure data; the harness decides what to run on them.

using UnderwaterAcoustics

# Canonical Munk profile [COA (5.104)], sampled for SampledField
function munk_ssp(; zmax=5000.0, n=101)
  z̃(z) = 2 * (z - 1300.0) / 1300.0
  c(z) = 1500.0 * (1.0 + 0.00737 * (z̃(z) - 1 + exp(-z̃(z))))
  zs = range(0.0, zmax; length=n)
  SampledField(c.(zs); z=-zs)
end

# 30-point Munk-like profile from test/test_raysolver.jl (seamount scenario)
function seamount_ssp()
  c = [1509.0, 1509.0, 1508.0, 1506.0, 1503.0, 1502.0, 1500.0, 1498.0, 1496.0,
       1495.0, 1492.0, 1491.0, 1490.0, 1489.0, 1488.0, 1487.0, 1487.0, 1487.0,
       1487.0, 1488.0, 1490.0, 1492.0, 1493.0, 1495.0, 1497.0, 1500.0, 1502.0,
       1505.0, 1508.0, 1511.0, 1515.0]
  z = 0.0:-100.0:-3000.0
  SampledField(c; z=z)
end

const SCENARIOS = [

  (name = "pekeris",
   description = "Isovelocity Pekeris waveguide: 20 m flat, sandy seabed, 1 kHz",
   tests = "straight-line rays, boundary reflection loss, coherent summation",
   env = UnderwaterEnvironment(bathymetry=20.0, soundspeed=1500.0, seabed=SandyClay),
   tx = AcousticSource(0.0, -5.0, 1000.0),
   grid = AcousticReceiverGrid2D(1.0:5.0:2000.0, -19.5:0.5:-0.5),
   rxs = [AcousticReceiver(500.0, -10.0), AcousticReceiver(1000.0, -5.0), AcousticReceiver(2000.0, -15.0)],
   referee = :pekeris),

  (name = "lloyd",
   description = "Lloyd mirror: 5 km deep isovelocity, near-surface source, 150 Hz",
   tests = "surface-image interference, geometric spreading, beam widths near free field",
   env = UnderwaterEnvironment(bathymetry=5000.0, soundspeed=1500.0,
                               seabed=FluidBoundary(water_density(), 1500.0, 0.0)),  # transparent bottom
   tx = AcousticSource(0.0, -10.0, 150.0),
   grid = AcousticReceiverGrid2D(10.0:10.0:3000.0, -200.0:2.0:-2.0),
   rxs = [AcousticReceiver(500.0, -50.0), AcousticReceiver(1500.0, -100.0), AcousticReceiver(3000.0, -30.0)],
   referee = :lloyd),

  (name = "munk50",
   description = "Munk profile deep water, 50 Hz, 100 km",
   tests = "refraction, convergence zones, caustics/KMAH, turning points",
   env = UnderwaterEnvironment(bathymetry=5000.0, soundspeed=munk_ssp(),
                               seabed=FluidBoundary(1600.0, 1600.0, dBperλ(0.5))),
   tx = AcousticSource(0.0, -1000.0, 50.0),
   grid = AcousticReceiverGrid2D(200.0:200.0:100_000.0, -5000.0:50.0:0.0),
   rxs = [AcousticReceiver(30_000.0, -1000.0), AcousticReceiver(60_000.0, -2000.0), AcousticReceiver(90_000.0, -800.0)],
   referee = :kraken),

  (name = "munk230",
   description = "Munk profile deep water, 230 Hz, 100 km",
   tests = "same as munk50 at higher frequency (finer interference structure)",
   env = UnderwaterEnvironment(bathymetry=5000.0, soundspeed=munk_ssp(),
                               seabed=FluidBoundary(1600.0, 1600.0, dBperλ(0.5))),
   tx = AcousticSource(0.0, -1000.0, 230.0),
   grid = AcousticReceiverGrid2D(200.0:200.0:100_000.0, -5000.0:50.0:0.0),
   rxs = [AcousticReceiver(30_000.0, -1000.0), AcousticReceiver(60_000.0, -2000.0), AcousticReceiver(90_000.0, -800.0)],
   referee = nothing),

  (name = "wedge",
   description = "Upslope wedge: 200 m → 20 m over 4 km, isovelocity, 100 Hz",
   tests = "range-dependent bathymetry, slope-corrected reflections, mode stripping",
   env = UnderwaterEnvironment(
     bathymetry=SampledField([200.0, 20.0, 20.0]; x=[0.0, 4000.0, 5000.0]),
     soundspeed=1500.0, seabed=SandyClay),
   tx = AcousticSource(0.0, -100.0, 100.0),
   grid = AcousticReceiverGrid2D(10.0:10.0:5000.0, -200.0:2.0:-2.0),
   rxs = [AcousticReceiver(1000.0, -50.0), AcousticReceiver(3000.0, -30.0), AcousticReceiver(4500.0, -10.0)],
   referee = nothing),

  (name = "seamount",
   description = "Seamount + Munk-like SSP, 230 Hz, 100 km (from test suite)",
   tests = "combined non-flat bathymetry + non-iso SSP, shadow zone behind seamount",
   env = UnderwaterEnvironment(
     bathymetry=SampledField([3000.0, 3000.0, 500.0, 3000.0, 3000.0];
                             x=[0.0, 10e3, 20e3, 30e3, 100e3]),
     soundspeed=seamount_ssp(),
     seabed=FluidBoundary(1600.0, 1600.0, dBperλ(0.5))),
   tx = AcousticSource(0.0, -18.0, 230.0),
   grid = AcousticReceiverGrid2D(200.0:200.0:100_000.0, -3000.0:30.0:0.0),
   rxs = [AcousticReceiver(15_000.0, -500.0), AcousticReceiver(40_000.0, -1000.0), AcousticReceiver(80_000.0, -500.0)],
   referee = nothing),

  (name = "downref",
   description = "Downward-refracting shallow gradient: 1500→1480 m/s over 40 m, silt, 500 Hz",
   tests = "many bottom bounces; reflection model and near-boundary beam accumulation",
   env = UnderwaterEnvironment(
     bathymetry=40.0,
     soundspeed=SampledField([1500.0, 1480.0]; z=[0.0, -40.0]),
     seabed=SandySilt),
   tx = AcousticSource(0.0, -10.0, 500.0),
   grid = AcousticReceiverGrid2D(5.0:5.0:3000.0, -39.0:1.0:-1.0),
   rxs = [AcousticReceiver(1000.0, -20.0), AcousticReceiver(2000.0, -35.0), AcousticReceiver(3000.0, -10.0)],
   referee = :kraken),

  (name = "duct",
   description = "Surface duct: upward-refracting near-surface gradient, 1 kHz",
   tests = "ducted propagation and leakage below the duct (ray-theory weak spot)",
   env = UnderwaterEnvironment(
     bathymetry=1000.0,
     soundspeed=SampledField([1500.0, 1505.0, 1490.0, 1485.0]; z=[0.0, -100.0, -400.0, -1000.0]),
     seabed=SandyClay),
   tx = AcousticSource(0.0, -30.0, 1000.0),
   grid = AcousticReceiverGrid2D(20.0:20.0:10_000.0, -500.0:5.0:-5.0),
   rxs = [AcousticReceiver(3000.0, -50.0), AcousticReceiver(6000.0, -30.0), AcousticReceiver(9000.0, -200.0)],
   referee = :kraken),

]

scenario(name::AbstractString) = only(filter(s -> s.name == name, SCENARIOS))
