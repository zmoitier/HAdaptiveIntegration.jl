# HAdaptiveIntegration

[![Lint workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Test workflow status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl)
[![Docs workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![code style: runic](https://img.shields.io/badge/code_style-%E1%9A%B1%E1%9A%A2%E1%9A%BE%E1%9B%81%E1%9A%B2-black)](https://github.com/fredrikekre/Runic.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev)

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
- Arbitrary-precision arithmetic, via a `ForwardDiff` extension that re-derives
  tabulated rules at higher precision (install `ForwardDiff` to enable it).
- Compatibility with [`Unitful.jl`](https://github.com/JuliaPhysics/Unitful.jl)
  quantities as integrand values and domain coordinates.

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
- `buffer` to reuse heap memory across repeated `integrate` calls (see
  `allocate_buffer`); useful for reducing allocations in a hot loop.
- `callback` to observe each estimate during refinement, receiving
  `(I, E, nb_subdiv, buffer)` on every step.

For full API details and advanced usage, see the
[stable documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable/) or the
[development documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev/).

## Related Packages and When to Use Them

`HAdaptiveIntegration` is inspired by
[`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl), which adaptively integrates
over orthotopes in any dimension. `HAdaptiveIntegration` adds support for simplices of
any dimension and uses tabulated cubatures for low-dimensional orthotopes, where it can
be competitive with `HCubature` using fewer function evaluations.

This package ships rules for arbitrary `d`-dimensional simplices and orthotopes, but
other packages may suit your case better:

- for `d=1` (where both `1`-simplex and `1`-orthotope reduce to a segment),
  [`QuadGK.jl`](https://github.com/JuliaMath/QuadGK.jl) is usually preferable;
- for medium-dimensional orthotopes,
  [`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl) may be faster;
- for large-dimensional simplices or orthotopes, deterministic adaptive cubature may
  become slow, so consider stochastic methods such as
  [`MCIntegration.jl`](https://github.com/numericalEFT/MCIntegration.jl) or
  [`Cuba.jl`](https://github.com/giordano/Cuba.jl).

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for how to
report a bug, ask a question, or submit a pull request.
