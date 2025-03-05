```@meta
CurrentModule = HAdaptiveIntegration
```

# HAdaptiveIntegration

Documentation for
[HAdaptiveIntegration](https://github.com/zmoitier/HAdaptiveIntegration.jl).

## Basic usage

The main function of this package is `integrate`, which computes the integral of a
[`Domain`](@ref)

```@example
using InteractiveUtils # hide
using HAdaptiveIntegration
HAdaptiveIntegration.LIST_DOMAIN
```

## Comparison with `HCubature.jl`

This package shares many similarities with `HCubature.jl`; there are, however, a few
important differences:

- `HAdaptiveIntegration.jl` uses tabulated embedded cubature rules such as the ones
  found in [Cubature.jl](https://www.google.com/?client=safari), whereas `HCubature.jl`
  implemented the Genz-Malik algorithm valid axis-aligned rectangles in any dimension.
- `HAdaptiveIntegration.jl` supports simplicies in low dimensions, whereas
  `HCubature.jl` supports axis-aligned rectangles in any dimension.

Let's start with a simple example using `HCubature`:

```@example hcubature-square
using HCubature, LinearAlgebra
a, b = (0.0, 0.0), (1.0,1.0)
const counter = Ref(0)
# f = x -> (counter[]+=1; 1 / (norm(x) + 1e-0))
f = x -> (counter[]+=1; cos(20*prod(x)))
I, E = hcubature(f, a, b)
println("I = $I, E = $E, counter = $(counter[])")
```

Now, let's do the same with `HAdaptiveIntegration`:

```@example hcubature-square
using HAdaptiveIntegration
domain = rectangle(a, b)
counter[] = 0
I, E = integrate(f, domain)
println("I = $I, E = $E, counter = $(counter[])")
```

Lets look at performance now:

```@example hcubature-square
using BenchmarkTools
counter[] = 0
b1 = @benchmark hcubature($f, $a, $b)
```

```@example hcubature-square
ec = HAdaptiveIntegration.embedded_cubature(HAdaptiveIntegration.SQUARE_CH21G25, Float64)
counter[] = 0
b2 = @benchmark integrate($f, $domain, $ec)
```

Let's do the same comparison for the 3d-cube.

```@example hcubature-cube
using HCubature, LinearAlgebra
a, b = (0.0, 0.0, 0.0), (1.0, 1.0, 1.0)
const counter = Ref(0)
f = x -> (counter[]+=1; cos(5*prod(x)))
I, E = hcubature(f, a, b)
println("I = $I, E = $E, counter = $(counter[])")
```

```@example hcubature-cube
using HAdaptiveIntegration
domain = cuboid(a, b)
counter[] = 0
I, E = integrate(f, domain)
println("I = $I, E = $E, counter = $(counter[])")
```

```@example hcubature-cube
using BenchmarkTools
counter[] = 0
b1 = @benchmark hcubature($f, $a, $b)
```

```@example hcubature-cube
ec = HAdaptiveIntegration.embedded_cubature(HAdaptiveIntegration.CUBE_BE65, Float64)
counter[] = 0
b2 = @benchmark integrate($f, $domain, $ec)
```

<!-- function Base.parse(T::Type{MultiFloat{Float64,N}}, str::String) where {N}
    return T(str)
end -->
