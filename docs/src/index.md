```@meta
CurrentModule = HAdaptiveIntegration
```

# HAdaptiveIntegration

## Overview

`HAdaptiveIntegration.jl` is a Julia package for approximating integrals of functions over
various predefined [`AbstractDomain`](@ref)s. It uses *embedded cubature* rules to build
error estimates, and refines the integration domain by splitting its mesh elements until a
certain tolerance is reached. Features include:

- Adaptive integration over **simplices and orthotope of any dimension**
- Use of **efficient (tabulated) cubatures** for low-dimensional simplices (triangle and
  tetrahedron) and orthotope (rectangle and cuboid)
- Support for custom cubature rules
- Arbitrary precision arithmetic

## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add HAdaptiveIntegration
```

Or, equivalently, via the Julia Pkg API

```julia
julia> import Pkg; Pkg.add("HAdaptiveIntegration")
```

## Basic usage

The main function exported by this package is [`integrate(f, Ω)`](@ref), which is used to
approximate

```math
I = \int_{\Omega} f(x) \, \mathrm{d}x
```

where ``\Omega \subset \mathbb{R}^d`` is a [`AbstractDomain`](@ref) object, and
``f \colon \mathbb{R}^d \to \mathbb{F}`` is a function. Here is a simple example: first, we
define a function,

```@example quickstart
using HAdaptiveIntegration

fct = x -> cis(sum(x)) / (sum(abs2, x) + 1e-2)
nothing # hide
```

!!! warning "Function signature"
    The function `f` must accept a single argument `x` which is an abstract vector of length
    `d`, the dimension of the domain (concretely, `f` is called through `f(::SVector)`). The
    return type `T` of `f` can be any type that supports the operations `+(T, T)`,
    `norm(T)`, and multiplication by a scalar (*e.g.* vectors or matrices).

`Domain`s are constructed using the following functions (see their respective docstrings for
more details):

- [`Triangle`](@ref): a triangle in 2D,
- [`Tetrahedron`](@ref): a tetrahedron in 3D,
- [`Simplex`](@ref): a simplex in arbitrary dimension,
- [`Rectangle`](@ref): a rectangle in 2D,
- [`Cuboid`](@ref): a cuboid in 3D,
- [`Orthotope`](@ref): an orthotope (hyperrectangle) in arbitrary dimension.

### Simplices

To integrate the above function over a ``d``-dimensional simplices (triangle, tetrahedron,
...), defined by their vertices, we can use

- Triangle

  ```@example quickstart
  I, E = integrate(fct, Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)))
  ```

  The result `I` is the integral of `f` over a triangle with vertices `(0,0)`, `(1,0)`, and
  `(0,1)`, and `E` is an error estimate.
- Tetrahedron

  ```@example quickstart
  I, E = integrate(fct, Tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
  ```

- ``4``-simplex

  ```@example quickstart
  I, E = integrate(
           fct,
           Simplex{Float64}((0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1));
           rtol=1e-4
         )
  ```

  The keyword arguments `atol` and `rtol` can be used to control the desired absolute and
  relative error tolerances, respectively.

### Orthotopes (*a.k.a.* hyperrectangle)

To integrate the same function over a ``d``-dimensional axis-aligned orthotope (rectangle,
cuboid, and hyperrectangle), defined by their low and high corners, we can use

- Rectangle

  ```@example quickstart
  I, E = integrate(fct, Rectangle((0.0, 0.0), (1.0, 1.0)))
  ```

- Cuboid

  ```@example quickstart
  I, E = integrate(fct, Cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
  ```

- ``4``-orthotope (hyperrectangle)

  ```@example quickstart
  I, E = integrate(fct, Orthotope((0.0, 0.0, 0.0, 0.0), (1.0, 1.0, 1.0, 1.0)); rtol=1e-4)
  ```

!!! tip "Related package to integration over an orthotope (hyperrectangle)"
    This package contains rule for an arbitrary ``d``-dimensional orthotope, however:
    - For ``d=1``, you may want to check
      [`QuadGk.jl`](https://github.com/JuliaMath/QuadGK.jl), as it is specialized to do
      adaptive integration over the segment.
    - For high ``d``, you may want to check
      [`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl), as it supports adaptive
      integration over arbitrarily high-dimensional axis-aligned orthotope.
    - For even larger ``d``, you may want to check
      [`MCIntegration.jl`](https://github.com/numericalEFT/MCIntegration.jl) or
      [`Cuba.jl`](https://github.com/giordano/Cuba.jl), as they use stochastic method.

## Going further

In the previous examples we covered the basic usage of the [`integrate`](@ref) function.
There are, however, other options that can be passed to `integrate` in order to customize
various aspects of the underlying algorithm (*e.g.* passing a different cubature rule, using
a buffer to avoid memory allocations, etc.). For more details, see the docstring of the
[`integrate`](@ref) function, as well as the next section on
[advanced usage](@ref advanced-usage).
