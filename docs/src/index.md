```@meta
CurrentModule = HAdaptiveIntegration
```

# HAdaptiveIntegration

## Overview

`HAdaptiveIntegration.jl` is a Julia package for approximating integrals of functions over
various predefined [`AbstractDomain`](@ref)s. It uses *embedded cubature* rules to build error
estimates, and refines the integration domain by splitting its mesh elements until a certain
tolerance is reached. Features include:

- Adaptive integration over **simplices of any dimension**
- Use of **efficient (tabulated) cubatures** for low-dimensional cuboids and simplices
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
``f \colon \mathbb{R}^d \to \mathbb{F}`` is a function. Here is a simple example:

```@example quickstart
using HAdaptiveIntegration
domain = triangle((0, 0), (1, 0), (0, 1))
fct    = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
I, E   = integrate(fct, domain)
```

The result `I` is the integral of `f` over a triangle with vertices `(0,0)`, `(1,0)`, and
`(0,1)`, and `E` is an error estimate.

!!! warning "Function signature"
    The function `f` must accept a single argument `x` which is an abstract vector of length
    `d`, the dimension of the domain (concretely, `f` is called through `f(::SVector)`). The
    return type `T` of `f` can be any type that supports the operations `+(T,T)`, `norm(T)`,
    and multiplication by a scalar (*e.g.* vectors or matrices).

The keyword arguments `atol` and `rtol` can be used to control the desired absolute and
relative error tolerances, respectively:

```@example quickstart
I, E = integrate(fct, domain; rtol = 1e-12)
```

Finally, to integrate the same function over a three-dimensional axis-aligned cuboid we can
use

```@example quickstart
domain = cuboid((0, 0, 0), (1, 1, 1))
I, E   = integrate(fct, domain)
```

`Domain`s are constructed using the following functions (see their respective docstrings for
more details):

- [`triangle`](@ref): a triangle in 2D,
- [`tetrahedron`](@ref): a tetrahedron in 3D,
- [`simplex`](@ref): a simplex in arbitrary dimension,
- [`rectangle`](@ref): a rectangle in 2D,
- [`cuboid`](@ref): a cuboid in 3D.

!!! tip "N-cuboids and `HCubature.jl`"
    If you are looking for a package that supports adaptive integration over arbitrarily
    high-dimensional axis-aligned cuboids, you may want to check out
    [`HCubature.jl`](https://github.com/JuliaMath/HCubature.jl).

## Going further

In the previous examples we covered the basic usage of the [`integrate`](@ref) function.
There are, however, other options that can be passed to `integrate` in order to customize
various aspects of the underlying algorithm (*e.g.* passing a different cubature rule, using
a buffer to avoid memory allocations, etc.). For more details, see the docstring of the
[`integrate`](@ref) function, as well as the next section on
[advanced usage](@ref advanced-usage).
