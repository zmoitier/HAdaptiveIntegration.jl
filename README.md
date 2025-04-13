# HAdaptiveIntegration

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev)
[![Test workflow status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

`HAdaptiveIntegration` is a Julia package designed for numerical integration over
multidimensional domains. It computes integrals of the form

```math
I = \int_{\Omega} f(x) \, \mathrm{d}x
```

where $f$ is any Julia function and $\Omega$ represents domains such as simplices and
orthotopes. The package employs an adaptive approach, dynamically refining the integration
domain as needed. It uses embedded quadrature rules to provide error estimates, aiming to
achieve high accuracy while minimizing function evaluations.

Features include:

- Adaptive integration over simplices and orthotope of **any dimension**,
- Utilization of **efficient tabulated cubatures** for low-dimensional simplices and
  orthotopes,
- Support for custom embedded cubature rules,
- Arbitrary precision arithmetic.

## Quick Examples

Here are simple examples of how to use `HAdaptiveIntegration` to compute an integral over
the supported shapes:

```julia
using HAdaptiveIntegration

# Define a function
f = x -> cis(sum(x)) / (sum(abs2, x) + 1e-2)

# Compute the integral and error estimate over a triangle and a rectangle
I, E = integrate(f, triangle((0, 0), (1, 0), (0, 1)))
I, E = integrate(f, rectangle((0, 0), (1, 1)))

# Compute the integral and error estimate over a tetrahedron and a cuboid
I, E = integrate(f, tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)))
I, E = integrate(f, cuboid((0, 0, 0), (1, 1, 1)))

```

The result `I` is the integral value and `E` the error estimate.

There are many options available for the `integrate` function, as well as other supported
integration domains. For more information, see the
[stable documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable/) or the
[latest development documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev/).

## Related packages

`HAdaptiveIntegration` draws inspiration from the
[`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl) package, which offers a similar
approach for integrating over orthotope in any dimension. The key differences are:

- `HAdaptiveIntegration` supports integration over simplices of any dimension, whereas
`HCubature` is focused on orthotopes.
- For low-dimensional orthotopes such as squares and cubes, `HAdaptiveIntegration` employs
  tabulated cubatures for enhanced efficiency. This allows it to achieve precision
  comparable to `HCubature` with fewer function evaluations for these domains.
