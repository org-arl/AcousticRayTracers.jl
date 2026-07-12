# RaySolver vs BELLHOP benchmarks

- `scenarios.jl` — scenario definitions (pure data; both models consume the same `UnderwaterEnvironment`).
- `harness.jl` — model-agnostic runners, metrics, convergence sweeps, timing, version-keyed caching.
- `run_benchmarks.jl` — populates `results/cache` (slow; run before rendering).
- `report.qmd` — the benchmark report; renders from the cache.

## Usage

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'   # first time
julia --project=benchmarks -t auto benchmarks/run_benchmarks.jl
quarto render benchmarks/report.qmd
```

Version pinning: `UnderwaterAcoustics#master` and `AcousticsToolbox#compat-uwa-0.8`
(RaySolver on this branch requires the UWA 0.8 ∂-field API; the registered
AcousticsToolbox caps UWA at 0.7). The results cache is keyed by the
AcousticRayTracers git SHA + UWA revision and auto-invalidates on change;
`FORCE_RERUN=1` forces recomputation.
