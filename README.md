[![doc-stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://org-arl.github.io/UnderwaterAcoustics.jl/raysolver.html)
[![CI](https://github.com/org-arl/AcousticRayTracers.jl/workflows/CI/badge.svg)](https://github.com/org-arl/AcousticRayTracers.jl/actions)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Codecov](https://codecov.io/gh/org-arl/AcousticRayTracers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/org-arl/AcousticRayTracers.jl)
[![ColPrac](https://img.shields.io/badge/ColPrac-contributing-blueviolet)](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md)

# AcousticRayTracers

`RaySolver` is a differentiable 2½D Gaussian beam tracer for use with [`UnderwaterAcoustics.jl`](https://github.com/org-arl/UnderwaterAcoustics.jl).
It is similar to [Bellhop](http://oalib.hlsresearch.com/AcousticsToolbox/), but fully written in Julia to be compatible with automatic differentiation (AD)
tools such as `ForwardDiff`.

`BellhopJL` is a native Julia port of the 2D Bellhop Gaussian beam / ray tracer
(following the [A-New-BellHope](https://github.com/A-New-BellHope/bellhop) bug-fixed
mirror of Bellhop 2022_4). It is a drop-in alternative to `AcousticsToolbox.Bellhop`
with no file I/O, and is likewise ForwardDiff-differentiable.

For information on how to use the models, see [documentation](https://org-arl.github.io/UnderwaterAcoustics.jl/raysolver.html).

## Licensing

This package is dual-licensed on a per-directory basis:

- Everything **except** `src/bellhopjl/` is licensed under the [MIT license](LICENSE).
- `src/bellhopjl/` (the `BellhopJL` solver) is a derivative work of the GPL-licensed
  Bellhop (© 1983–2022 Michael B. Porter; A-New-BellHope changes © 2021–2023 The
  Regents of the University of California) and is licensed under
  [GPL-3.0-or-later](LICENSE-GPL-3.0). Each file in that directory carries an
  `SPDX-License-Identifier: GPL-3.0-or-later` header.

**Note:** because the combined package includes the GPL-3.0 `BellhopJL` component,
distributing the package as a whole is subject to the terms of the GPL-3.0. If you
need MIT-only terms, strip `src/bellhopjl/` (the rest of the package does not depend
on it).

## Contributing

Contributions in the form of bug reports, feature requests, ideas/suggestions, bug fixes, code enhancements, and documentation updates are most welcome. Please read [contribution guidelines](https://github.com/org-arl/UnderwaterAcoustics.jl/blob/master/CONTRIBUTING.md) if you wish to start contributing.

The scopes active in this repository are:
- **raysolv**: RaySolver
- **bellhopjl**: BellhopJL
