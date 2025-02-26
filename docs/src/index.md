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
subtypes(HAdaptiveIntegration.Domain)
```

## Comparisson with `HCubature.jl`

This package shares many similarities with `HCubature.jl`; there are, however, a few
important differences:

- `HAdaptiveIntegration.jl` uses tabulated embedded cubature rules such as the ones
  found in [Cubature.jl](https://www.google.com/?client=safari), whereas `HCubature.jl`
  implemented the Genz-Malik algorithm valid axis-aligned rectangles in any dimension.
- `HAdaptiveIntegration.jl` supports simplicies in low dimensions, whereas
  `HCubature.jl` supports axis-aligned rectangles in any dimension.

Let's start with a simple example using `HCubature`:

```@example hcubature
using HCubature, LinearAlgebra
a, b = (0.0, 0.0), (1.0,1.0)
const counter = Ref(0)
# f = x -> (counter[]+=1; 1 / (norm(x) + 1e-0))
f = x -> (counter[]+=1; 1+cos(20*prod(x)))
I, E = hcubature(f, a, b)
println("I = $I, E = $E, counter = $(counter[])")
```

Now, let's do the same with `HAdaptiveIntegration`:

```@example hcubature
import HAdaptiveIntegration as HAI
domain = HAI.rectangle(a, b)
counter[] = 0
I, E = HAI.integrate(f, domain)
println("I = $I, E = $E, counter = $(counter[])")
```

Lets look at performance now:

```@example hcubature
using BenchmarkTools
counter[] = 0
b1 = @benchmark hcubature($f, $a, $b)
```

```@example hcubature
counter[] = 0
b2 = @benchmark HAI.integrate($f, $domain)
```
