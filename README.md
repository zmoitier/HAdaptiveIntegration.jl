# HAdaptiveIntegration

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev)
[![Test workflow status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/zmoitier/HAdaptiveIntegration.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/zmoitier/HAdaptiveIntegration.jl)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

`HAdaptiveIntegration` is a Julia package designed for numerical integration over multi-dimensional domains. It computes integrals of the form

```math
I = \int_{\Omega} f(x) \, dx
```

where $f$ is any Julia function and $\Omega$ represents domains such as simplices and cuboids. The package employs an adaptive approach, dynamically refining the integration domain as needed. It uses embedded quadrature rules to provide error estimates, aiming to achieve high accuracy while minimizing function evaluations.

## Features

- Adaptive integration over **simplices of any dimension**
- Utilization of **efficient tabulated cubatures** for low-dimensional cuboids and simplices
- Support for custom cubature rules
- Arbitrary precision arithmetic

## Quick Example

```julia
using HAdaptiveIntegration

# Define a triangle domain and function to integrate
tri = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)

# Compute the integral and error estimate over the triangle
I, E = integrate(f, tri)
```

The result `I` is the integral value, and `E` the error estimate.

There are many options available for the `integrate` function, as well as other supported
integration domains. For more information, see the [stable
documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/stable/) or the [latest
development documentation](https://zmoitier.github.io/HAdaptiveIntegration.jl/dev/).

## Related packages

`HAdaptiveIntegration` draws inspiration from the
[`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl) package, which offers a similar
approach for integrating over hyperrectangles in any dimension. The key differences are:

- `HAdaptiveIntegration` supports integration over simplices of any dimension, whereas
  `HCubature` is focused on hyperrectangles.
- For low-dimensional domains such as squares, cubes, and triangles, `HAdaptiveIntegration`
  employs tabulated cubatures for enhanced efficiency. This allows it to achieve precision
  comparable to `HCubature` with fewer function evaluations for these domains.
