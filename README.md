# HAdaptiveIntegration

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev)
[![Test workflow status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

`HAdaptiveIntegration` is a Julia package for adaptive numerical integration on
multidimensional simplices and orthotopes. It computes integrals of the form

```math
I = \int_{\Omega} f(x) \, \mathrm{d}x
```

where $f$ is a Julia function and $\Omega$ is the integration domain. The algorithm
adaptively subdivides the domain and uses embedded cubature rules to estimate error,
targeting high accuracy with fewer function evaluations.

**Key Features:**
- Adaptive integration over simplices and orthotopes of arbitrary dimension.
- Efficient tabulated cubature rules for low-dimensional simplices and orthotopes.
- Support for custom embedded cubature rules.
- Arbitrary-precision arithmetic.

## Installation

```julia
using Pkg
Pkg.add("HAdaptiveIntegration")
```

## Quick Start

```julia
using HAdaptiveIntegration

f = x -> cis(sum(x)) / (sum(abs2, x) + 1e-2)

# Segment
I, E = integrate(f, Segment(0, 1))

# Triangle and rectangle
I, E = integrate(f, Triangle((0, 0), (1, 0), (0, 1)))
I, E = integrate(f, Rectangle((0, 0), (1, 1)))

# Tetrahedron and cuboid
I, E = integrate(f, Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
I, E = integrate(f, Cuboid((0, 0, 0), (1, 1, 1)))

# 4-simplex and 4-orthotope
I, E = integrate(
    f,
    Simplex((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1));
    rtol=1e-4,
)
I, E = integrate(f, Orthotope((0, 0, 0, 0), (1, 1, 1, 1)); rtol=1e-4)
```

`I` is the integral estimate and `E` is the error estimate.

## Practical `integrate` Options

Common keyword arguments:

- `atol` (absolute tolerance) and `rtol` (relative tolerance) to control stopping
  tolerances.
- `maxsubdiv` to cap the number of refinements.

For full API details and advanced usage, see the
[stable documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable/) or the
[latest development documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev/).

## Related Packages and When to Use Them

`HAdaptiveIntegration` draws inspiration from the
[`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl) package, which offers a similar
approach for integrating over orthotopes in any dimension. Key differences:

- `HAdaptiveIntegration` supports integration over simplices of any dimension, whereas
  `HCubature` is focused on orthotopes.
- For low-dimensional orthotopes such as squares and cubes, `HAdaptiveIntegration` employs
  tabulated cubatures for efficiency. In these domains it can be competitive with
  `HCubature` while using fewer function evaluations.

This package includes rules for arbitrary `d`-dimensional simplices and orthotopes, but:

- for `d=1` (where both `1`-simplex and `1`-orthotope reduce to a segment),
  [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) is usually preferable;
- for medium-dimensional orthotopes,
  [`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl) may be faster;
- for large-dimensional simplices or orthotopes, deterministic adaptive cubature may become
  slow, so it may be better to use stochastic methods such as
  [`MCIntegration.jl`](https://github.com/numericalEFT/MCIntegration.jl) or
  [`Cuba.jl`](https://github.com/giordano/Cuba.jl).
