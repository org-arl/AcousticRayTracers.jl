# BELLHOP internals, compared to RaySolver

Audience: the RaySolver (AcousticRayTracers.jl) maintainer.

Sources (all local, read-only):

- **BELLHOP Fortran** — A-New-BellHope mirror of the 2022_4 release:
  `/Users/mandar/Projects/BellhopJL.jl/references/bellhop` (2D program is
  `Bellhop/bellhop.f90` + modules; `misc/` holds shared Acoustics-Toolbox code).
  Comments prefixed `LP:` in that source are annotations/fixes by the
  A-New-BellHope maintainers, and its `README.md` documents the numerical
  reproducibility problems of the original OALIB code.
- **AcousticsToolbox.jl wrapper** —
  `/Users/mandar/Projects/BellhopJL.jl/references/AcousticsToolbox.jl/src/{bellhop.jl,common.jl}`.
- **RaySolver** —
  `/Users/mandar/Projects/AcousticRayTracers-benchmark/src/RaySolver.jl` (~730
  lines in the benchmark checkout, including backscatter/scatterer support).
- **BellhopJL.jl** native port + `PORTING_NOTES.md`
  (`/Users/mandar/Projects/BellhopJL.jl`) — used for verified constants and
  sign conventions.

Throughout, `file:subroutine` refers to the BELLHOP source tree and
`RaySolver.jl:NN` to the benchmark copy of RaySolver.

Recommendations are graded **ADOPT / PARTIAL / DON'T ADOPT**, always subject to
RaySolver's two hard constraints: (a) the whole trace must remain
ForwardDiff-differentiable (no hard-coded `Float64`, no branchy state machines
that break dual propagation), and (b) the trace is a SciML `ODEProblem` with
`VectorContinuousCallback` event handling, not a hand-rolled stepper.

---

## 1. Top-level architecture

### 1.1 Program flow (`bellhop.f90:BELLHOP` / `BellhopCore`)

BELLHOP is a batch program keyed off a *file root*: it reads `root.env` plus
optional companion files, and writes one of `root.shd` (TL), `root.arr`
(arrivals), `root.ray` (ray/eigenray trace), always echoing to `root.prt`.
The driver (`bellhop.f90:1–148`) does, in order:

1. `ReadEnvironmentBell.f90:ReadEnvironment` — parse the `.env` file and
   initialize the SSP (§2, §8).
2. `bdryMod.f90:ReadATI` / `ReadBTY` — optional altimetry/bathymetry files
   (`root.ati`, `root.bty`), triggered by option letters `TopOpt(5:5)` and
   `BotOpt(2:2)` being `'~'`/`'*'` (`ReadEnvironmentBell.f90:334, 98`).
3. `misc/RefCoef.f90:ReadReflectionCoefficient` — optional tabulated top/bottom
   reflection coefficients (`root.trc`/`root.brc`), triggered by boundary
   condition letter `'F'`.
4. `misc/beampattern.f90:ReadPat` — optional source beam pattern (`root.sbp`),
   triggered by `RunType(3:3) == '*'`.
5. `BellhopCore` (`bellhop.f90:152–358`): for each source depth, zero the
   pressure/arrival matrix, then loop over the launch-angle fan
   (`DeclinationAngle` loop, `bellhop.f90:267–322`). For each angle: look up the
   source beam pattern amplitude `Amp0` by linear interpolation
   (`bellhop.f90:272–279`), apply the Lloyd-mirror factor for semicoherent runs
   (`bellhop.f90:282–283`, §7), call `TraceRay2D`, then dispatch on
   `Beam%Type(1:1)` to one of six influence functions (`bellhop.f90:300–318`,
   §5). Finally `ScalePressure` and write to disk.

Key structural contrast with RaySolver: BELLHOP **always traces the full fan
once and splats each beam onto the receiver grid**; eigenrays and arrivals are
just different "influence" backends over the same trace. RaySolver has two
distinct modes: an *eigenray search* (`arrivals`, `RaySolver.jl:83–134`) that
root-finds on the receiver-depth error over launch angle, and a *Gaussian-beam
field* (`acoustic_field(::AcousticReceiverGrid2D)`, `RaySolver.jl:151–243`)
that is architecturally the same as BELLHOP's influence pass. BELLHOP has no
root-finding anywhere — its "eigenrays" (`RunType 'E'`) are simply every ray
whose beam window contains the receiver, written out as ray files
(`influence.f90:ApplyContribution`, case `'E'`).

**Recommendation — DON'T ADOPT the monolithic structure.** RaySolver's split is
a feature: exact eigenrays (needed for differentiable impulse responses) are
something BELLHOP cannot produce at all. Keep it.

### 1.2 Default fan size and step size

- `Beam%deltas == 0` ⇒ `deltas = depth/10` (`bellhop.f90:169–173`). RaySolver's
  `ds` default is identical: `minimum(bathymetry)/10` (`RaySolver.jl:94, 162`).
- Auto beam count (`angleMod.f90:38–52`):
  `Nalpha = max(⌊0.3·Rmax·f/1500⌋, 300)` and at least
  `π/atan2(depth, 10·Rmax)`; 50 for ray-trace plots. Separately BellhopCore
  *warns* if `Nalpha < NalphaOpt = 2 + (αmax−αmin)/√(c/(6·f·Rmax))`
  (`bellhop.f90:257–263`). Note both auto-counts are **frequency- and
  range-dependent**.
- The AcousticsToolbox.jl wrapper overrides this for arrivals runs with a flat
  0.05° spacing: `nbeams = (max_angle−min_angle)/0.05° + 1` ≈ **3201 beams**
  for the default ±80° fan (`bellhop.jl:53`).
- RaySolver's auto heuristics are geometric only: arrivals use
  `⌈16·(Δα)/atan(h,R)⌉` (`RaySolver.jl:92`), the beam field uses
  `clamp(⌈8·Δα·R/Δz⌉, 100, 1000)` (`RaySolver.jl:160`) — no frequency term,
  and typically one to two orders of magnitude fewer beams (benchmark finding
  (iv), §6).

**Recommendation — PARTIAL.** RaySolver's eigenray search converges each
arrival exactly, so it does not need BELLHOP's dense fan for *accuracy of an
individual arrival*; but the fan density controls *which* arrivals are found
(bracketing) and the smoothness of the beam field. Two cheap adoptions:
(1) add the frequency-aware `NalphaOpt`-style floor
(`2 + Δα·√(6·f·Rmax/c)`) to the field heuristic, since Gaussian-beam fringe
resolution genuinely scales with `√f·R`; (2) for arrivals, raise the auto count
(or document the risk) in ducted/multi-turning-point profiles where 16 rays per
`atan(h,R)` under-brackets closely spaced eigenray families. Do not adopt the
wrapper's 0.05° default wholesale — with root-finding it wastes ~10× work.

---

## 2. SSP interpolation (`sspMod.f90`)

`sspMod.f90:EvaluateSSP` dispatches on `SSP%Type` = `TopOpt(1:1)` and returns
`c, cimag, gradc(2), crr, crz, czz, rho` — i.e. the value **and all second
derivatives** the dynamic-ray equations need, plus the imaginary sound speed
that encodes volume attenuation (§7). Options (`ReadEnvironmentBell.f90:258–285`):

| Letter | Routine | Interpolant | `czz` supplied |
|---|---|---|---|
| `C` | `sspMod.f90:cLinear` | piecewise-linear in c | 0 within a layer |
| `N` | `sspMod.f90:n2Linear` | piecewise-linear in n²=1/c² | `czz = 3·cz²/c` (`sspMod.f90:216`) |
| `S` | `sspMod.f90:cCubic` | cubic spline, **not-a-knot** ends (`iBCBeg=iBCEnd=0` → `misc/splinec.f90:CSPLINE`) | analytic |
| `P` | `sspMod.f90:cPCHIP` | monotone piecewise-cubic Hermite (`misc/pchipMod.f90`) | analytic (`sspMod.f90:301`) |
| `Q` | `sspMod.f90:Quad` | bilinear in (r,z) from `root.ssp`; `crr = czz = 0`, `crz` finite-differenced across range columns | 0 |
| `H` | `sspMod.f90:Hexahedral` | trilinear 3D (BELLHOP3D only) | 0 |
| `A` | `sspMod.f90:Analytic` | hard-coded Munk-like profile | analytic |

Details that matter:

- Coefficients are **precomputed at read time** (`Task='INI'`,
  `sspMod.f90:ReadSSP`); per-step evaluation is a segment lookup plus a cubic
  evaluation. Segment lookup (`sspMod.f90:UpdateDepthSegmentT`, 867–892) is
  stateful (`iSegz` module variable) and *tangent-direction-aware*: a ray
  sitting exactly on a knot is assigned to the segment it is about to enter.
  This hysteresis is what makes the exact knot-landing of the stepper (§3)
  well-defined.
- `SSP%c` is complex (`sspMod.f90:32`): `misc/AttenMod.f90:CRCI` converts the
  user attenuation to an imaginary part at read time; `cimag` then feeds the
  `τ` update `w/complex(c, cimag)` in the stepper (`Step.f90:93`), so volume
  attenuation is integrated *along the exact ray path* as `Im(ωτ)`.
- Density is interpolated linearly in every profile routine but is only used
  in halfspace reflection (water-side ρ is taken as the env value).
- **Piecewise-linear layers make ray arcs exact**: within a `C`-linear layer,
  `czz = 0` and the ray equations have circular-arc solutions; BELLHOP's
  layer-clamped RK2 (§3) commits a step exactly at each knot and then applies
  a discrete p-jump (§3.2), which together handle the distributional `∂²c`
  (delta at the knot) *exactly*.

### RaySolver equivalent

RaySolver takes `c(z)`, `∂c`, `∂²c` from the environment via
`UnderwaterAcoustics.∂` (`RaySolver.jl:392–394`) — the interpolation model is
whatever the env's `soundspeed` field provides (constant, `SampledFieldZ` with
linear/cubic interpolation, or an arbitrary differentiable function). This is
strictly more general in one direction (any AD-able `c(z)`, including neural
SSPs) and narrower in another: it is range-independent by construction
(`RaySolver.jl:344–351`, `_∂u` sets `dξ/ds = 0`), while BELLHOP's `Q` option
and `crr`/`crz` terms handle range-dependent SSPs.

The Gaussian-beam `q, p` ODE in `_∂u` uses `c̄ = ∂²c(z)·ξ²` as `cnn/c²`
(`RaySolver.jl:349`), the range-independent specialization of BELLHOP's
`cnn0_csq0 = crr·t₂² − 2·crz·t₁·t₂ + czz·t₁²` (`Step.f90:50`; note the
misleading Fortran name — it is *not* divided by c², see PORTING_NOTES).

**Recommendations.**
- **N²-linear (`N`): DON'T ADOPT.** It exists for closed-form legacy reasons;
  a c-linear or spline SSP through UnderwaterAcoustics covers all use cases,
  and adding another interpolant class belongs in UnderwaterAcoustics, not the
  tracer.
- **PCHIP (`P`): PARTIAL — via UnderwaterAcoustics.** Monotone interpolation
  avoids the spurious ducts a not-a-knot spline can invent between coarse CTD
  points; that's a real modeling win. But it should be a `SampledFieldZ`
  interpolation mode (with AD-able coefficients) consumed transparently by
  RaySolver's existing `∂` interface, not RaySolver code. Note even the
  AcousticsToolbox.jl wrapper never selects PCHIP (`common.jl:54–56` maps only
  to `'C'`/`'S'`/`'Q'`).
- **Range-dependent SSP (`Q`): PARTIAL, long-term.** Supporting it means
  generalizing `_∂u` to carry `∂c/∂r` (state derivative on ξ becomes nonzero)
  and `cnn` in full quadratic form. The SciML formulation handles this
  naturally and differentiably; the cost is API (env must expose `c(r,z)`
  derivatives) and re-verifying `∂raysolver`. Worth doing only when a use case
  arrives; the `_∂u` comment already flags the assumption.
- **Complex-c volume attenuation along the path: ADOPT the idea** (see §7).

---

## 3. The stepper (`Step.f90`)

### 3.1 `Step2D` — layer-clamped RK2

`Step.f90:Step2D` (11–135) is a modified midpoint/Heun second-order step of the
state `(x, t, p, q, τ)` where `t` is the 1/c-scaled tangent:

1. **Phase 1** (Euler half step): evaluate SSP at `ray0`, set `h = deltas`,
   call `ReduceStep2D` to clamp `h` (below), take a half step `halfh = h/2`
   for position, tangent, p, q (`Step.f90:54–64`).
2. **Phase 2**: re-evaluate SSP at the midpoint, clamp again with the midpoint
   tangent, then blend derivatives with weights `w1 = h/(2·halfh)`, `w0=1−w1`
   (`Step.f90:86–93`) — this preserves 2nd-order accuracy when the clamp
   shrank `h` between phases.
3. **Commit** via `StepToBdry2D` (`Step.f90:98`), which recomputes the minimal
   clamp for the blended tangent and, crucially, **snaps the landing position
   exactly onto the boundary value stored in memory** (e.g.
   `x2(2) = SSP%z(iSegz0)`, `Step.f90:248`) and returns `topRefl`/`botRefl`
   flags. (This snap + flags is the A-New-BellHope overhaul; the original OALIB
   code landed "randomly" on either side of every boundary — see README.md
   "Boundary stepping changes" — which is why every step of every OALIB ray is
   an edge case.)

`ReduceStep2D` (`Step.f90:139–225`) clamps `h` to the first of: SSP-layer
crossing in depth (`hInt`), beam-box r/z exit (`hBoxr`, `hBoxz`), top/bottom
boundary crossing via signed distance along the outward normal (`hTop`,
`hBot`), and boundary-*segment* range crossing (`hSeg`, including SSP range
columns when `Type=='Q'`). A floor `h ≥ 1e-6·deltas` forces motion and counts
sequential small steps (`iSmallStepCtr`, not used for termination in 2D).

### 3.2 The SSP-interface p-jump

After committing a step, if the layer index changed (`Step.f90:117–133`),
BELLHOP applies the discrete curvature jump

```
RM = +t₁/t₂  (depth crossing)  or  −t₂/t₁  (range crossing)
RN = RM·(2·cnjump − RM·csjump)/c₂ ;  p ← p − q·RN
```

with `cnjump/csjump` = jumps of `∇c` along the ray normal/tangent. This is the
distributional `∂²c` (delta function at a c-linear knot) that a smooth ODE
integrator cannot see. **This is exactly the mechanism behind benchmark finding
(i)**: RaySolver's `_∂u` integrates `∂²c` pointwise, so with a piecewise-linear
SSP it missed the delta contribution to `p` at every knot — fixed by
AcousticRayTracers **PR #26** (SSP-knot transmission correction), which adds
the equivalent jump at knot crossings.

### RaySolver equivalent and recommendations

RaySolver replaces the whole stepper with `solve(ODEProblem, Tsit5/Rosenbrock23)`
plus a 6-component `VectorContinuousCallback` (`_check_ray!`,
`RaySolver.jl:357–364`) for surface/bottom/max-range/turn-back/scatterer/
receiver-crossing events, with adaptive error control (`solver_tol`) instead of
a fixed `deltas`.

- **Fixed-h RK2: DON'T ADOPT.** Adaptive high-order SciML solvers dominate on
  accuracy-per-cost and are the differentiable path. BELLHOP's stepper exists
  because it predates good ODE infrastructure and because exact knot landing
  is its interface-handling strategy.
- **Layer clamping / exact knot events: PARTIAL (largely done).** The
  ODE-callback analogue of `ReduceStep2D`'s `hInt` is an event function on SSP
  knot depths; with PR #26's p-jump this reproduces BELLHOP's exact per-segment
  treatment while staying differentiable. If PR #26 currently detects knots via
  tolerance rather than a callback, consider promoting knot crossings to a
  proper `ContinuousCallback` — it removes solver-tolerance sensitivity at
  knots the same way `StepToBdry2D`'s snap removed edge-case randomness in
  Fortran.
- **Small-step floor / `iSmallStepCtr`: DON'T ADOPT.** It is a symptom of the
  fixed stepper; SciML's `dtmin` plus RaySolver's `atol` clamp of the ray
  position at reflection (`RaySolver.jl:536`) already prevent stalls, and the
  StepsizeLimiter anti-tunneling guard (`RaySolver.jl:432`) handles the
  scatterer case.
- **Beam box: PARTIAL.** BELLHOP terminates rays on a user `Box%r/Box%z`
  (wrapper writes `1.01×` grid extents, `common.jl:164`). RaySolver's
  `rmax`/`rmin` cover range; there is no depth box, but bathymetry/altimetry
  events make one redundant in 2D. No action needed.

---

## 4. Ray trace loop and reflection

### 4.1 `TraceRay2D` (`bellhop.f90:437–609`)

Initial conditions: `t = (cosα, sinα)/c`, `p = (1,0)`, `q = (0,1)`
(`bellhop.f90:464–470`); geometric runs zero `q(2)` since only the real ray
tube is tracked. The loop: `Step2D` → update altimetry/bathymetry segment
(tangent-aware, `bellhop.f90:519–543`) → if `topRefl`/`botRefl`, interpolate
the boundary normal (constant per segment for `'L'` piecewise-linear
boundaries; node-interpolated for the `'C'` curvilinear option,
`bellhop.f90:556–563`) and call `Reflect2D`. Termination
(`bellhop.f90:592–603`): outside the beam box, `Amp < 0.005`, two successive
points outside a boundary, or step-storage exhausted (`MaxN = 100000`).
Reflection stores **two coincident ray points** (incoming and outgoing
tangent), which the influence functions later skip as duplicates.

### 4.2 `Reflect2D` (`bellhop.f90:633–805`)

1. Specular tangent update `t ← t − 2(t·n̂)n̂` (`bellhop.f90:659`).
2. **Dynamic-ray curvature correction** (Müller 1984): boundary-curvature term
   `RN = 2κ/c²/Th` plus the soundspeed-gradient jump term
   `RM(2·cnjump − RM·csjump)/c²`, sign-flipped for the top boundary, applied as
   `p ← p + q·RN` (`bellhop.f90:661–698`). `Beam%Type(3:3)` can double (`'D'`)
   or zero (`'Z'`) the correction; default `'S'` = single.
3. **Reflection coefficient**, dispatch on `HS%BC` (`bellhop.f90:702–803`):
   - `'R'` rigid: unchanged. `'V'` vacuum: phase += π.
   - `'F'` file: interpolate the tabulated |R|, φ table
     (`misc/RefCoef.f90:InterpolateReflectionCoefficient` — linear in angle,
     **R set to 0 outside the tabulated domain**, `RefCoef.f90:146–158`).
   - `'A'`/`'G'` halfspace: analytic fluid/elastic R from `kzP`, `kzS`
     (`bellhop.f90:715–742`), with a compiler-branch fix on `SQRT` of negative
     reals (`bellhop.f90:736–737`). Rays with `|R| < 1e-5` are killed.
   - Optional **beam displacement** (Tindle shift, `Beam%Type(4:4)=='S'`,
     "Seongil's version", `bellhop.f90:763–797`): translates the reflection
     point along the bottom, adds `pdelta` to τ and a `q` width change. The
     code itself notes it is unreliable with segment-crossing logic.

Amplitude/phase are carried as separate real `Amp`, `Phase` on the ray
(`bellhopMod` ray2DPt), i.e. |R| and arg R accumulate separately.

### RaySolver equivalent and recommendations

RaySolver's reflection handling is in `_trace` (`RaySolver.jl:556–581`):
callback event 1/2 fires on the boundary, `reflection_coef(...)` from
UnderwaterAcoustics supplies the **exact** complex R (accumulated into the
complex phasor `A` — cleaner than BELLHOP's split Amp/Phase), the angle is
slope-corrected via `θ ← −θ ± 2·atan(boundary′)`, and `_reflect_pq`
(`RaySolver.jl:451–463`) applies the same Müller curvature correction with
true boundary curvature `κ = f″/(1+f′²)^{3/2}` obtained from the field
derivative API (`∂(..., :x, :x)`) — this is the curved-boundary work shipped in
AcousticRayTracers PR #25. Note the sign convention difference: RaySolver's
`RN = −2κ/(c²Th) + ...` with κ positive convex-into-water and n̂ pointing
*into* the water, vs BELLHOP's `+2κ/c²/Th` with outward normals and a top-side
sign flip — verified equivalent (see `_reflect_pq` docstring and
PORTING_NOTES "Reflect2D").

- **Curvature correction: already adopted** (PR #25). One divergence:
  RaySolver skips the correction for grazing incidence `|Th|·c ≤ 1e-3`
  (`RaySolver.jl:453`) where BELLHOP lets it blow up (or the user selects
  `'Z'`). Keep RaySolver's guard — it is also the AD-friendly choice.
- **Reflection-coefficient physics: DON'T ADOPT BELLHOP's.** UnderwaterAcoustics'
  `reflection_coef` is exact and differentiable; BELLHOP's is equivalent for
  the fluid halfspace but is fed *translated* geoacoustic parameters when
  driven via the wrapper — that translation (density in g/cm³, attenuation
  re-derived through `CRCI` in dB/λ) is the source of benchmark finding (ii):
  −3 to −4 dB discrepancies on 2–3-bounce wedge rays near the intromission
  angle come from the wrapper's env translation, not from RaySolver.
- **Tabulated R (TRC/BRC): PARTIAL.** Useful when a measured bottom-loss curve
  is all you have. The right home is a UnderwaterAcoustics boundary type
  (`TabulatedBoundary`?) whose `reflection_coef` interpolates (AD-able if the
  table values can carry duals); RaySolver would then get it for free. Low
  effort, medium value.
- **Beam displacement: DON'T ADOPT (now).** Even BELLHOP's own comments flag it
  as unreliable; it matters mainly for low-frequency shallow-water TL fine
  structure. Revisit only if a benchmark shows a need.
- **`|R| < 1e-5` ray kill and `Amp < 0.005` termination: already covered** by
  `min_amplitude` (`RaySolver.jl:586–589`) — but see §6(v): RaySolver's
  `min_amplitude = 1e-6` on `|A|/|q|` culls weak multipath that BELLHOP still
  bins into arrivals. Consider making the arrivals-path cull threshold
  independent of (and lower than) the field-path one, or exposing it more
  prominently, if weak-multipath completeness matters to users.

---

## 5. Influence functions (`influence.f90`) and normalization

BELLHOP separates "trace the beam" from "deposit the beam". Six backends,
selected by `Beam%Type(1:1)` (`bellhop.f90:300–318`):

### 5.1 Geometric hat, Cartesian (`influence.f90:InfluenceGeoHatCart`, 420–522) — the default

Per ray segment `[iS−1, iS]`: build the unit tangent/normal from the segment
chord, skip duplicate points, track caustic phase (`IncPhaseIfCaustic`,
`influence.f90:900–909`: +π/2 whenever `q` crosses zero, sign convention
`qleq0=true`). Walk a **single receiver-range pointer** monotonically (bracket
`[rA, rB)`, bump toward `rB`, `influence.f90:505–517` — O(#steps + #receivers)
per ray, direction-aware so backward-travelling rays work). For each bracketed
receiver: project onto the segment (`s`, `n`), interpolate
`q = q(iS−1) + s·dq`, beam half-width `RadiusMax = |q/q0|` with
`q0 = c₀/Δα` (`influence.f90:430`), and if `n < RadiusMax` deposit

```
Amp = Ratio1·√(c/|q|)·rayAmp · (RadiusMax − n)/RadiusMax     (hat taper)
field += Amp · exp(−i(ω·delay − phase))
```

`Ratio1 = √|cosα|` for a point source (`influence.f90:442`). A vertical-window
pre-cull (`zmin/zmax` from segment endpoints ± RadiusMax) is applied only for
shallow rays (`|rayt(1)| > 0.5`).

### 5.2 Geometric Gaussian, Cartesian (`InfluenceGeoGaussianCart`, 526–651) — `Beam%Type 'B'`

Same walk, but width `σ = |q/q0|` floored by
`max(σ, min(0.2·f·Re τ, π·λ))` (`influence.f90:581, 616`) — a near-field/
minimum-width clamp RaySolver lacks — window `n < 4σ`, taper
`W = exp(−½(n/σ)²)/(σ·|q0/q|)`, and `Ratio1` gains a `1/√(2π)` so a fan of
Gaussians sums to the free-space field (`influence.f90:552–556`). The
shallow-angle threshold is `0.50001` (A-New-BellHope moved it off the exact
60° edge, README "Shallow angle range").

### 5.3 Ray-centered variants (`InfluenceGeoHatRayCen`, 310–416)

Instead of interpolating along the ray at fixed receiver ranges, they solve for
the intersection of each receiver depth line with the *ray-normal* coordinate
(`nB = (zR − z(iS))/znV(iS)`), which handles steep/vertical rays better.
Caustic phase is precomputed per step (`KMAHphase`, `influence.f90:325–333`)
— the A-New-BellHope fix for receiver-placement-dependent phase (README
"Precomputing phase in 2D hat ray-centered influence").

### 5.4 Cerveny beams (`InfluenceCervenyRayCen` 19–167, `InfluenceCervenyCart` 171–306) — brief

The original physics-based Gaussian beams: complex
`qVB = q₁ + ε·q₂`, `Γ = p/q`, field
`∝ Amp·√(c·|ε|/q)·exp(−i(ω(τ + ½Γn²) − phase))`, with ε chosen by
`bellhop.f90:PickEpsilon` (space-filling `ε = i·½ω·(2/(ω/c·Δα))²`,
minimum-width, or WKB; several codepaths leave ε undefined — flagged UB,
`bellhop.f90:376`). Supports image beams in the boundary mirrors
(`Beam%Nimage` ≤ 3) and a KMAH branch-cut tracker on complex q
(`influence.f90:BranchCut`, 775–799). TL-only (no arrivals/eigenrays).

### 5.5 SGB (`InfluenceSGB`, 685–771)

Bucker's simple Gaussian beams; contains an acknowledged bug (assumes every
step has length `deltas`, `influence.f90:729–733`). Legacy; ignore.

### 5.6 `ApplyContribution` and `ScalePressure`

`ApplyContribution` (`influence.f90:655–681`) branches the same interpolated
(Amp, delay, phase) into: eigenray file write, `AddArr`, coherent sum, or
incoherent/semicoherent power sum (`(const·e^{Im ωτ})²·W`, with an extra
`√(2π)` for Gaussian). `ScalePressure` (`influence.f90:803–841`) applies the
final normalization: `const = −1` for geometric beams (the sign the wrapper
cancels by negating shd pressures, `common.jl:312` and PORTING_NOTES
"ScalePressure"), `−Δα·√f/c` for Cerveny; `√(Re U)` converts incoherent
intensity to pressure; point sources get `1/√r` cylindrical spreading (0 at
r = 0), line sources `−4√π·const`.

### RaySolver equivalent and recommendations

RaySolver's `beam!` closure (`RaySolver.jl:178–220`) is a Gaussian
(hat-less) Cartesian influence: per segment, find the receiver-grid subset
within `4·max(|q1|,|q2|)·δθ` (`_findall_rxgrid2d` — O(1) index arithmetic on
the regular grid, better than BELLHOP's pointer walk for grids), compute the
CPA projection `α`, width `W = |q·δθ|` [COA (3.74)], amplitude
`A = C1·√(δθ·|cosθ₀|·cₛ/(c₀·r·W))` [COA (3.76)], taper `exp(−(n/W)²)`, phase
`cispi(2f·t)` with KMAH already folded into the phasor at trace time
(`RaySolver.jl:533`). The COA-(3.76) normalization plays the role of
`ScalePressure` (there is no separate pass; `1/√r` is inside the per-deposit
formula, and the `√π` Gaussian-sum factor is in `C1`, `RaySolver.jl:190`).

- **Hat beams: DON'T ADOPT.** Hats are BELLHOP's default largely for speed and
  historical continuity; Gaussians are smoother, and RaySolver's tests/users
  are calibrated to them. (If exact BELLHOP mimicry is wanted, BellhopJL.jl
  already exists.)
- **Gaussian minimum-width clamp: ADOPT.** BELLHOP's
  `σ ← max(σ, min(0.2·f·Re τ, π·λ))` prevents unphysically narrow beams near
  the source and near caustics where `q → 0` (where RaySolver's `W = |q·δθ|`
  collapses and the `1/√W` amplitude blows up; the KMAH phasor handles the
  phase but not the amplitude singularity). Direct, differentiable, ~3 lines
  in `beam!`. This is the single highest-value influence-function adoption.
- **Ray-centered influence: PARTIAL/LOW.** RaySolver's segment-CPA projection
  already handles steep segments (it projects onto the segment, not onto
  vertical lines), so the main ray-centered advantage is moot. Skip unless
  vertical-ray artifacts show up in benchmarks.
- **Cerveny beams: DON'T ADOPT.** More physics, but complex-q tracing doubles
  the ODE state, PickEpsilon is heuristic (and partly UB in the reference),
  TL-only, and the geometric-beam formulation is the community default. Not
  worth the AD and maintenance surface.
- **Below-boundary deposits — benchmark finding (iii).** BELLHOP's influence
  functions deposit into *any* receiver inside the beam window, including
  receivers below the seafloor (the halfspace is a legitimate part of its
  domain; rays are killed but beams still overlap sub-bottom grid points).
  RaySolver clips deposits to the water column (rays never leave it, and the
  grid subset is centered on in-water segments). **DON'T ADOPT** BELLHOP's
  behavior — sub-bottom "field" from a water-column beam sum is not physical —
  but benchmark comparisons must mask sub-bottom receivers before diffing.
- **Line sources (`RunType(4:4)=='X'`): DON'T ADOPT** unless a user asks;
  it's a normalization-only change (`Ratio1 = 1`, `−4√π` factor) and would be
  easy to add later behind a solver flag.

---

## 6. Caustics/KMAH, arrivals, eigenrays

### 6.1 KMAH tracking

BELLHOP (geometric beams): phase accumulator += π/2 on every sign change of
real `q` along the ray (`influence.f90:IncPhaseIfCaustic/IsAtCaustic`,
900–927), plus a per-deposit correction: if the *receiver-interpolated* q has
crossed relative to the segment start, `FinalPhase` (`influence.f90:874–896`)
adds π/2 — and (probable bug, kept for parity in BellhopJL) **discards the
accumulated ray phase** in that case; Gaussian reads the ray phase from the
current point, hat from the previous (another kept inconsistency). Cerveny
beams use the complex branch-cut tracker instead (`BranchCut`).

RaySolver: counts sign changes of `u[i][7]` (q) across saved points per
segment and multiplies the phasor by `cis(−π/2·kmah)` (`RaySolver.jl:506–510,
533`), also per-crossing for backscatter signatures (`RaySolver.jl:517–529`).

**Recommendation — PARTIAL.** RaySolver's approach is cleaner and avoids the
Fortran inconsistencies, but it detects crossings **at saved points**, so with
large `ds` (or `save_everystep=false` in the `ds ≤ 0` path,
`RaySolver.jl:436`) a caustic between saved points where q dips through zero
and back could be missed, and the deposit-time interpolated q in `beam!` does
not get BELLHOP's `FinalPhase`-style sub-segment correction. Two options:
(a) add q = 0 as a (non-terminal) continuous callback — exact, differentiable,
cheap; (b) add the sign-of-interpolated-q check in `beam!`. Option (a) is the
principled fix; skip BELLHOP's phase-discard quirk.

### 6.2 Arrivals accumulation (`ArrMod.f90:AddArr`, 23–87) — benchmark finding (v)

Every beam deposit at a receiver becomes an arrival: `(A, Phase, complex
delay, src/rcvr angles, bounce counts)`. Merging: if `ω·|Δdelay| < 0.05` and
`|Δphase| < 0.05` versus the **last stored arrival only** (flagged bug, kept),
amplitude-weight-merge delay and angles, sum amplitudes; else append; when the
per-receiver store (`MaxNArr = 20000000/(NRz·NRr)`, min 10,
`bellhop.f90:219`) is full, replace the weakest arrival if the new one is
stronger. `WriteArrivalsASCII` applies the `1/√r` spreading factor at write
time (`ArrMod.f90:105–113`). Because *every* beam within a window contributes,
BELLHOP reports many weak micropath/bracketing arrivals; adjacent-beam pairs
hitting the same receiver merge into one. RaySolver instead reports one exact
eigenray per root of Δz(θ), and its `min_amplitude=1e-6` cull inside `_trace`
(`RaySolver.jl:589`) additionally terminates weak rays — hence finding (v):
BELLHOP's arrival binning reports weak multipath that RaySolver drops.

**Recommendation — DON'T ADOPT binning; PARTIAL on the cull.** Binned arrivals
have amplitude-weighted (i.e. slightly wrong) delays and beam-count-dependent
counts; exact eigenrays are strictly better for impulse responses and for AD.
But the completeness gap is real: expose/loosen the arrivals-path amplitude
cull (it is a *trace termination* threshold, so it also truncates paths that
would later brighten — rare, but possible after a caustic), and document that
`min_amplitude` trades weak-multipath completeness for speed.

### 6.3 Eigenrays

BELLHOP `RunType 'E'`: any ray whose beam window covers the receiver gets its
whole trajectory written to the ray file (`ApplyContribution` case `'E'`).
The wrapper then *matches* these approximate eigenrays to arrivals by bounce
counts and nearest launch angle to attach paths (`bellhop.jl:58–68`) — an
inherently fuzzy join. RaySolver's bracketing + `IntervalNonlinearProblem` /
turning-point `NonlinearProblem` search (`RaySolver.jl:107–132`) is exact.
**DON'T ADOPT.**

### 6.4 Cross-referenced benchmark findings (given, from the study in this directory)

1. **SSP-knot transmission correction** — BELLHOP handles piecewise-linear
   SSPs exactly per segment (layer-clamped steps + `Step.f90:117–133` p-jump);
   RaySolver's smooth ODE missed the `∂²c` delta at knots. Fixed by
   AcousticRayTracers **PR #26**. (§2, §3.2.)
2. **Intromission-angle reflection differences** — the wedge-scenario −3 to
   −4 dB deltas on 2–3-bounce rays trace to the AcousticsToolbox wrapper's
   geoacoustic env translation (g/cm³ density, dB/λ attenuation via `CRCI`)
   versus UnderwaterAcoustics' exact `reflection_coef`; not a RaySolver defect.
   (§4.2.)
3. **Sub-seafloor receivers** — BELLHOP deposits field into below-seafloor
   (halfspace) grid points; RaySolver clips to the water column. Mask before
   comparing. (§5.)
4. **Beam counts** — BELLHOP-via-wrapper defaults to 0.05° spacing (~3201
   beams) for arrivals; RaySolver's auto heuristic is far sparser
   (`RaySolver.jl:92, 160`). Sparse fans are fine for root-found eigenrays but
   can under-bracket and under-resolve fields. (§1.2.)
5. **Weak multipath** — BELLHOP's `AddArr` binning reports arrivals that
   RaySolver's `min_amplitude=1e-6` culling drops. (§6.2.)

---

## 7. Attenuation and run modes

### 7.1 `misc/AttenMod.f90:CRCI` (21–123)

Converts (speed, attenuation) to complex `c` at env-read time. Units
`AttenUnit(1:1)`: `N` Np/m, `F` dB/(m·kHz), `M` dB/m, `W` dB/λ
(`αT = α·f/(8.6858896·c)` Np/m), `Q` quality factor, `L` loss parameter; then
`cimag = αT·c²/ω`. `AttenUnit(2:2)` adds volume attenuation: `T` Thorp
(`AttenMod.f90:94–97`, the JKPS Eq. 1.34 variant in dB/km → Np/m), `F`
Francois–Garrison (`AttenMod.f90:98–99, Franc_Garr` 127–176; T, S, pH, z̄ read
from the env — `ReadEnvironmentBell.f90:313`), `B` biological layers. Because
`cimag` rides in the SSP, volume attenuation accumulates as `Im(ωτ)` along the
exact path through `Step.f90:93`.

RaySolver applies **nominal** absorption once per arrival:
`absorption(f, D, salinity, mid_temp, hmax/2)` on total path length D at
mid-depth (`RaySolver.jl:594`), and per-segment `exp(logγ·s̄)` in the beam
field (`RaySolver.jl:171–186`).

**Recommendation — PARTIAL.** For a 2½D water column the nominal
D-times-α(f, mid-depth) approximation errs only through the depth/temperature
dependence of α, usually ≪ 0.1 dB/km of error; not worth complexifying the ODE
state. If depth-dependent absorption ever matters (strong thermoclines at
high f), the differentiable route is adding `∫α(z(s))ds` as an 8th ODE state —
straightforward — rather than complex c.

### 7.2 Coherent / incoherent / semicoherent

BELLHOP `RunType(1:1)`: `C` coherent phasor sum; `I` incoherent
(`|contri|²` summed, `√(Re U)` at the end, `influence.f90:673–678, 825`);
`S` semicoherent = incoherent **plus a Lloyd-mirror source factor**
`Amp0 ·= √2·|sin(ω·zs·sinα/c)|` applied at launch (`bellhop.f90:282–283`) —
i.e. incoherent TL with coherent surface-image interference at the source.
RaySolver: `:coherent`/`:incoherent` in both field paths (`RaySolver.jl:139–146,
191, 213–215, 241`); no semicoherent.

**Recommendation — ADOPT semicoherent (cheap).** It is literally one
`√2·|sin(2πf·zs·sinθ₀/c₀)|` factor on `C1` in `beam!` under a
`mode === :semicoherent` branch of the existing incoherent path; trivially
differentiable; and the wrapper exposes it (`bellhop.jl:78`), so model
intercomparisons in the Julia ecosystem already exercise it.

---

## 8. ENV file format reference (as BELLHOP reads it, and what the wrapper writes)

Read by `ReadEnvironmentBell.f90:ReadEnvironment` (list-directed Fortran READs,
so `/` terminates a line and repeats defaults). Line by line, with the exact
output of `AcousticsToolbox.jl/src/common.jl:_write_env` noted:

1. **Title** (quoted string). Wrapper: temp-dir suffix (`common.jl:41`).
2. **Frequency** in Hz (`ReadEnvironmentBell.f90:52`). Wrapper: `%0.6f`
   (`common.jl:43`).
3. **NMedia** — must be 1 for BELLHOP (`ReadEnvironmentBell.f90:57`).
   Wrapper: 1 (multilayer counts are for Kraken; `common.jl:44–45`).
4. **TopOpt**, up to 6 letters (`ReadTopOpt`, 239–350):
   - (1:1) SSP interpolant: `N|C|P|S|Q|H|A` (§2). Wrapper: `'C'` default,
     `'S'` iff `SampledFieldZ` with `CubicSpline()`, `'Q'` iff
     `SampledFieldXZ` (`common.jl:54–56`) — **never `N` or `P`**.
   - (2:2) top BC: `V|R|A|G|F|W|P`. Wrapper: `R`/`V`/`A` from the surface
     boundary (`common.jl:57`) — **never `F`** (no TRC).
   - (3:4) attenuation unit + volume attenuation: wrapper always writes
     `"WF"` = dB/λ + **Francois–Garrison** (`common.jl:58`), so the next env
     line must be `T S pH z̄` — wrapper writes
     `temperature salinity pH waterdepth/2` (`common.jl:67`).
   - (5:5) `~`/`*` = read `root.ati`. Wrapper: `*` iff altimetry is
     range-dependent (`common.jl:60–63`).
   - (6:6) `I` = development options. Wrapper: never.
5. **Top halfspace line** (only if top BC `'A'`): `z cp cs rho alphaI betaI`.
   Wrapper: `0.0 c ρ/1000 0.0 in_dBperλ(δ) /` pattern (`common.jl:68–70`).
6. **`NPts Sigma Depth`** (`ReadEnvironmentBell.f90:73`) — NPts and Sigma are
   unused by BELLHOP; Depth is the bottom depth and tells `ReadSSP` when to
   stop. Wrapper: `0 0.0 waterdepth` (`common.jl:76`).
7. **SSP block**: lines of `z alphaR betaR rho alphaI betaI`
   (`sspMod.f90:ReadSSP`, 811–865), terminated when `z == Depth`; `/` lets the
   wrapper write only `z c /` per line. Wrapper (`common.jl:78–100`): constant
   SSP → two nodes at 0 and waterdepth; `SampledFieldZ` → the sample grid with
   a `−waterdepth` node appended and an early break at waterdepth;
   `SampledFieldXZ` → z-column at r=0 plus a companion `model.ssp` written by
   `_create_ssp_file` (`common.jl:201–216`; §8-SSP below); other callables →
   resampled at `_recommend_len` (λ/2 spacing, clamped 25–1000,
   `common.jl:178–182`).
8. **BotOpt + Sigma**: `'A|R|V|G|F' [~|*] sigma`
   (`ReadEnvironmentBell.f90:93–107`). Wrapper: `R`/`V`/`A` from the seabed
   (elastic collapsed to fluid when `cₛ=0`), `*` appended iff bathymetry is
   range-dependent (writes `model.bty`), then σ (`common.jl:114–127`) —
   **never `F` or `G`**.
9. **Bottom halfspace line** (BC `'A'`): `depth cp cs rho/1000 δp[dB/λ] δs`
   (`common.jl:128–134`).
10. **Source depths**: `NSz` then values `/` (`ReadSzRz`). Wrapper: 1 depth,
    `−z(tx)` (`common.jl:139`).
11. **Receiver depths**: `NRz` then values `/`. Wrapper: `−zrange` reversed to
    ascending, recording `zrev` (`common.jl:141–154`).
12. **Receiver ranges** in **km**: `NRr` then values `/`. Wrapper: `xrange/1000`
    ascending, recording `xrev`.
13. **RunType**, up to 7 letters (`ReadRunType`, 354–436):
    `(1) R|E|I|S|C|A|a` task; `(2) C|R|S|b|B|g|G` beam type;
    `(3) *` = read SBP file; `(4) R` point / `X` line source;
    `(5) R` rectilinear / `I` irregular grid; `(6) 2|3` dimensionality;
    `(7) S` beam shift → `Beam%Type(4:4)`
    (`ReadEnvironmentBell.f90:159`). Wrapper writes
    `'<task><bcode>[S]'` with bcode `G` (geometric hat) or `B` (Gaussian)
    only — `:cartesian`/`:ray_centered` Cerveny codes exist in the code but
    are not constructible (constructor restricts to
    `(:geometric, :gaussian)`, `bellhop.jl:20`) — and `S` iff `beam_shift`
    (`common.jl:160–161`). Task codes used: `C|I|S` (field), `A` (arrivals),
    `E` (eigenray paths), `R` (`rays()` debug API) (`bellhop.jl:54–119`).
    **Never `'*'` (SBP), never `X`, `I`, `3`.**
14. **Nalpha** then **angle pair** `α1 α2 /` in degrees, *positive down*
    (`angleMod.f90:ReadRayElevationAngles`; `Nalpha 0` = auto; a negative
    3rd entry `-999.9` triggers `SubTab` linear fill). Wrapper:
    `max(nbeams, pm.nbeams)` and `−max_angle −min_angle` (sign flip for the
    up-positive UA convention; `common.jl:162–163`).
15. **`deltas Box%z Box%r`** — step (m, 0 = auto depth/10), box depth (m), box
    range (**km**) (`ReadEnvironmentBell.f90:147–155`). Wrapper:
    `0.0 1.01·waterdepth 1.01·maxr/1000` (`common.jl:164`).
16. **Cerveny-only extra lines** (`Beam%Type(1:1) ∈ {R,C}`):
    `'<width><curv>' epsMult rLoop` and `Nimage iBeamWindow Component`
    (`ReadEnvironmentBell.f90:196–229`). Never written by the wrapper.

### Companion files

- **BTY / ATI** (`bdryMod.f90:ReadBTY/ReadATI`): first line interpolation type
  `'L'` (piecewise linear) or `'C'` (curvilinear — node-averaged normals and
  `kappa = dφ/ds` curvature, `bdryMod.f90:295–332`); optionally `'LL'` with
  per-segment geoacoustics (`btyType(2:2)=='L'`, cp/cs/rho per segment,
  converted via CRCI in `bellhop.f90:181–197`). Then `N` and `r[km] z[m]`
  pairs; boundaries are extended to ±√huge/1e5 beyond the data
  (PORTING_NOTES "Boundary normals"). Wrapper `_create_alt_bathy_file`
  (`common.jl:184–199`): writes `'L'` for `SampledFieldX` with `Linear()`
  interp, `'C'` otherwise, or resamples a functional bathymetry at λ/2 onto a
  `'L'` grid; altimetry values are negated (`(q,p) -> -value(q,p)`,
  `common.jl:62`). Never writes geoacoustic (`'LL'`) segments.
- **SSP (range-dependent, Type `'Q'`)** (`sspMod.f90:Quad`, 368–488):
  `Nr`, the profile ranges in km, then `NPts×Nr` sound speeds row per depth.
  Depth sampling must match the env-file z nodes. Wrapper:
  `_create_ssp_file` (`common.jl:201–216`).
- **TRC/BRC** (`misc/RefCoef.f90`): `N` then `theta[deg] |R| phi[deg]` triples,
  monotone in theta; also `.irc` precalculated internal RC for BC `'P'`.
  Not written by the wrapper.
- **SBP** (`misc/beampattern.f90`): `N` then `angle[deg] level[dB]` pairs,
  linearly interpolated at launch (`bellhop.f90:272–279`). Not written by the
  wrapper.

---

## 9. Feature-gap matrices

### 9.1 BELLHOP features not in RaySolver

| Feature (BELLHOP ref) | Value for RaySolver | Effort | Verdict |
|---|---|---|---|
| Nx2D / 3D (`bellhop3D.f90`, `Step3DMod`, `influence3D`) | High for 3D bathymetry problems, but a different product | Very high; 3D dynamic rays + AD | Don't adopt (separate package if ever) |
| N²-linear SSP (`sspMod.f90:n2Linear`) | Low | Low (UWA-side) | Don't adopt |
| PCHIP SSP (`sspMod.f90:cPCHIP`) | Medium (no spurious ducts) | Low–medium, belongs in UnderwaterAcoustics | Partial (upstream) |
| Range-dependent SSP `'Q'` (`sspMod.f90:Quad`) | Medium–high | Medium (`_∂u` generalization + env API) | Partial, long-term |
| Curvilinear boundary interp `'C'` (`bdryMod.f90:295`) | Low — RaySolver already gets smooth κ from `∂(bathymetry,:x,:x)` (PR #25), which is *better* | — | Already superseded |
| Per-segment boundary geoacoustics `'LL'` (`bdryMod.f90:86, bellhop.f90:484–496`) | Medium (range-dependent seabeds) | Medium; env API question | Partial (needs UWA range-dependent boundary type) |
| Tabulated reflection coefficients TRC/BRC (`RefCoef.f90`) | Medium (measured bottom loss) | Low as a UWA boundary type | Partial (upstream) |
| Source beam patterns SBP (`beampattern.f90`, `bellhop.f90:272–279`) | Medium (directional sources) | Low: multiply launch amplitude by an interpolated pattern; AD-able | Adopt if UWA grows a directional-source API |
| Semicoherent mode + Lloyd mirror (`bellhop.f90:282–283`) | Medium | Very low (§7.2) | **Adopt** |
| Gaussian minimum-width clamp (`influence.f90:581`) | High (caustic amplitude blowup) | Very low | **Adopt** (§5) |
| Hat beams / ray-centered influence (`influence.f90:310–522`) | Low | Medium | Don't adopt |
| Cerveny beams + PickEpsilon (`influence.f90:19–306`) | Low | High, AD-hostile | Don't adopt |
| Beam displacement on reflection (`bellhop.f90:763–797`) | Low (self-flagged unreliable) | Medium | Don't adopt |
| Line sources (`RunType(4:4)=='X'`, `ScalePressure`) | Low | Low | Defer |
| Irregular receiver grid `'I'` | Low (`arrivals` per-receiver covers it) | — | Don't adopt |
| Complex-c path-integrated volume attenuation (`AttenMod.f90:CRCI`, `Step.f90:93`) | Low–medium | Low (extra ODE state) if ever needed | Defer (§7.1) |
| Arrival binning/merging (`ArrMod.f90:AddArr`) | Negative (fuzzier than root-found eigenrays) | — | Don't adopt |
| Ray-file / eigenray-file outputs (`WriteRay.f90`) | N/A (RaySolver returns paths in-memory) | — | Superseded |
| Biological attenuation layers (`AttenUnit(2:2)=='B'`) | Low | Low–medium | Defer |
| Grain-size seabed `'G'` (`ReadEnvironmentBell.f90:493–539`) | Low (UWA users specify geoacoustics directly) | Low as a UWA helper | Don't adopt |
| Multiple source depths per run (`SourceDepth` loop) | Low (call the solver per tx) | — | Don't adopt |

### 9.2 BELLHOP features not exposed by the AcousticsToolbox.jl wrapper

Verified against `common.jl:_write_env` / `bellhop.jl` (line refs above):

| Feature | Wrapper status | Value / effort to expose |
|---|---|---|
| ATI altimetry files | **Exposed** (written iff altimetry range-dependent, `common.jl:60–63`) — listed for the record | — |
| TRC/BRC tabulated reflection coefficients | Never written (BC letters limited to `R|V|A`, `common.jl:57, 122`) | Medium value / low effort (write `.trc/.brc` + BC `'F'`) |
| Source beam pattern SBP (`RunType(3:3)='*'`) | Never written (`common.jl:161` writes only task+bcode+shift) | Medium / low |
| SSP interpolation choice `N`, `P` | Not selectable (`common.jl:54–56` maps only C/S/Q) | Low / trivial (a kwarg) |
| Cerveny beam types `C`/`R` (+ eps/Nimage lines) | Dead code: bcodes exist in `common.jl:160` but constructor forbids them (`bellhop.jl:20`), and the extra env lines are never written (would crash BELLHOP) | Low value / medium effort |
| SGB beams `'S'`, hat ray-centered `'g'` | Not selectable | Low / trivial |
| Ray-file output for arbitrary receivers (`RunType 'R'`) | Partially: `rays(pm, tx, max_range)` debug API only (`bellhop.jl:115–121`) | — |
| Semicoherent mode | **Exposed** (`mode=:semicoherent` → task `'S'`, `bellhop.jl:78`) | — |
| Line source `'X'` | Never (RunType position 4 left default `'R'`) | Low / trivial |
| Irregular grid `'I'`, 3D `'3'` | Never | Low / N-A |
| Beam box override | Not exposed; hard-coded `1.01×` extents (`common.jl:164`) | Low / trivial |
| `deltas` step-size override | Not exposed; always 0 = auto (`common.jl:164`) | Medium (accuracy studies) / trivial |
| Single-beam trace (`iSingle_alpha`, `TopOpt(6:6)='I'`) | Never | Debug-only |
| Attenuation unit / volume-attenuation choice | Hard-coded `"WF"` (dB/λ + Francois–Garrison, `common.jl:58`) — notably **not Thorp** | Low / trivial |
| Grain-size seabed `'G'`, biological attenuation `'B'` | Never | Low |
| Per-segment bathymetry geoacoustics (`'LL'` BTY) | Never (`_create_alt_bathy_file` writes only `'L'`/`'C'`, `common.jl:186–193`) | Medium / medium |
| Multiple source depths | Never (one tx per run) | Low |

---

## 10. Summary of recommendations

**Adopt now (cheap, high value):**
1. Gaussian beam minimum-width clamp `σ ← max(σ, min(0.2·f·Re τ, π·λ))` in
   `beam!` (§5.2) — fixes caustic amplitude blowup.
2. `:semicoherent` mode (Lloyd-mirror factor on the incoherent path) (§7.2).
3. q = 0 continuous callback for exact, save-grid-independent KMAH detection
   (§6.1).
4. Frequency-aware floor on the auto beam count (§1.2).

**Adopt via UnderwaterAcoustics (upstream, not RaySolver code):** PCHIP SSP
interpolation; tabulated-reflection-coefficient boundary type; (eventually)
directional sources / beam patterns and range-dependent seabed properties.

**Already adopted / superseded:** curved-boundary p-q reflection correction
(PR #25 ≙ `Reflect2D` Müller terms, with a better κ source than BELLHOP's
curvilinear-BTY option); SSP-knot p-jump (PR #26 ≙ `Step.f90:117–133`);
amplitude-based ray kill (`min_amplitude` ≙ `Amp < 0.005` / `|R| < 1e-5`).

**Don't adopt:** fixed-step layer-clamped RK2; hat/Cerveny/SGB influence;
arrival binning (`AddArr`); beam displacement; window-based "eigenrays";
below-seafloor field deposits; N²-linear SSP; 3D.

**Benchmark-comparison hygiene:** mask sub-seafloor receivers (finding iii);
match beam counts or expect fringe-level coherent-TL differences (finding iv,
and the A-New-BellHope README notes 601-vs-600-beam sensitivity of median
0.22 dB / p95 2.6 dB); expect BELLHOP weak-multipath arrivals below RaySolver's
cull (finding v); attribute intromission-angle reflection deltas to the
wrapper's env translation, not the tracer (finding ii).
