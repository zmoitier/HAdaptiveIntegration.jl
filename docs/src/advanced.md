```@meta
CurrentModule = HAdaptiveIntegration
```

# [Advanced usage](@id advanced-usage)

We now cover the options available for the `integrate` function.

## Buffering

When calling `integrate(f, domain)`, the package allocates memory for storing the various
subregions that are generated during the adaptive integration process. Here is what it looks
like in practice:

```@example buffering
using HAdaptiveIntegration
using BenchmarkTools
t = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
@benchmark integrate($f, $t)
```

While the overhead associated to these (small) allocations are usually negligible, there are
circumstances where one may want to avoid allocations altogether. This can achieved by
passing a buffer to the `integrate` using [`allocate_buffer`](@ref):

```@example buffering
using HAdaptiveIntegration: allocate_buffer, integrate, triangle
buffer = allocate_buffer(f, t)
integrate(f,t; buffer)
@benchmark integrate($f, $t; buffer = $buffer)
```

## Embedded cubature formulas

By default, when calling `integrate(f, domain)`, the package uses a default embedded
cubature formula for the given `domain` by calling [`default_embedded_cubature`](@ref).
Although these are generally good choices, you can also specify a custom embedded
cubature formula by passing it as the third argument to `integrate`. For example, in the
case of a triangle, the package defaults to a Radon-Laurie embedded cubature formula of
orders 5 and 8 (see the `rules_simplex.jl` file for more details). If you want e.g. to use an
embedded cubature based on the [`GrundmannMoeller`](@ref) rule of order 13, you can do

```@example embedded-cubature
using HAdaptiveIntegration: GrundmannMoeller, embedded_cubature, integrate, triangle
t = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
ec = embedded_cubature(GrundmannMoeller(2,13), Float64)
I,E = integrate(f, t; embedded_cubature = ec)
```

!!! tip "Available embedded cubature formulas"
    You can find a list of available embedded cubature formulas in the `rule_simplex.jl`
    file and `rule_orthotope.jl` files.

To add a custom embedded quadrature for a given domain, you must write a constructor e.g.
`my_custom_cubature(args...)` that returns a valid [`EmbeddedCubature`](@ref) object (see
`embedded_cubature.jl` for some examples on how this is done). PRs with new schemes are more
than welcome!

## Subdivision strategies

The package uses a default subdivision strategy for the given `domain` by calling
[`default_subdivision`](@ref). For example, by default triangles are subidivided
into 4 smaller triangles by connecting the midpoints of the edges:

```@example default-subdivision
using HAdaptiveIntegration
t = triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
subdiv_algo = HAdaptiveIntegration.default_subdivision(t)
```

Here are the subdivided triangles:

```@example default-subdivision
subdiv_algo(t)
```

But is is also possible (and maybe desirable) to split the triangle into 2 smaller triangles
instead. The following function accomplishes this:

```@example default-subdivision
using HAdaptiveIntegration: subdivide_triangle2
subdivide_triangle2(t)
```

Passing `subdivide_triangle2` as the `subdiv_algo` to `integrate` will use this instead of
the default:

```@example default-subdivision
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
I, E = integrate(f, t; subdiv_algo = subdivide_triangle2)
```

Which subdivision strategy is best depends on the function being integrated; for the example
presented above, it turns out the default strategy is more efficient!
