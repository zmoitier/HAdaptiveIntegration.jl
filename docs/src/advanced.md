```@meta
CurrentModule = HAdaptiveIntegration
```

# [Advanced usage](@id advanced-usage)

We now cover the options available for the [`integrate`](@ref) function.

## [Reduce memory allocations](@id reduce-mem-alloc)

When calling `integrate(f, domain)`, the package allocates memory for storing the various subregions that are generated during the adaptive integration process. Here is what it looks like in practice:

```@example buffering
using HAdaptiveIntegration
using BenchmarkTools

t = Triangle((0, 0), (1, 0), (0, 1))
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
@benchmark integrate($f, $t)
```

While the overhead associated with these (small) allocations is usually negligible, there are circumstances where one may want to avoid allocations altogether. This can be achieved by passing a buffer to the [`integrate`](@ref) using [`allocate_buffer`](@ref):

```@example buffering
using HAdaptiveIntegration: allocate_buffer

buffer = allocate_buffer(f, t)
integrate(f, t; buffer)
b = @benchmark integrate($f, $t; buffer = $buffer)
@assert b.allocs == 0  # hide
b  # hide
```

Provided evaluating `f` does not allocate, and the `buffer` has a sufficiently large capacity, `integrate` will not allocate memory during the integration process, as shown in the benchmark above.

!!! note "When to use a buffer"
    Buffer pre-allocation is useful in **hot loops** where `integrate` is called repeatedly, *e.g.*, in optimization or parameter fitting. For one-off integrations, the overhead is negligible.

## [Track convergence progress](@id callback-fct)

The `callback` keyword argument allows you to monitor the progress of the adaptive integration. The callback function is called after each subdivision and receives:

| Argument    | Type      | Description                               |
| ----------- | --------- | ----------------------------------------- |
| `I`         | `Float64` | Current integral estimate                 |
| `E`         | `Float64` | Error estimate (absolute, not relative)   |
| `nb_subdiv` | `Int`     | Number of subdivisions performed so far   |
| `buffer`    |           | Internal buffer (passed for advanced use) |

!!! note
    The integration stops when `E ≤ max(atol, rtol * |I|)`. Use `atol` and `rtol` keyword arguments to control the tolerance.

Here is a practical example to print the history convergence:

```@example callback
using HAdaptiveIntegration

t = Triangle((0, 0), (1, 0), (0, 1))
f = x -> 1 / (x[1]^2 + x[2]^2 + 0.5)

history = @NamedTuple{I::Float64, E::Float64, nb_subdiv::Int}[]

I, E = integrate(f, t; callback = (I, E, nb_subdiv, _) -> push!(history, (; I, E, nb_subdiv)))

using Printf
@printf("  %4s | %16s | %12s\n", "step", "I", "E")
@printf("  %4s-+-%16s-+-%12s\n", "----", "----------------", "------------")
foreach(history) do h
    @printf("  %4d | %16.12f | %12.4e\n", h.nb_subdiv, h.I, h.E)
end
```

## Choose custom cubature rules

By default, when calling `integrate(f, domain)`, the package uses a default embedded cubature formula for the given `domain` by calling [`default_rule`](@ref). Although these are generally good choices, you can also specify a custom embedded cubature formula by passing it as a keyword argument to `integrate`. For example, in the case of a triangle, the package defaults to a Radon-Laurie embedded cubature formula of high order 8 and low order 5 (see [`RadonLaurie`](@ref)). If you want *e.g.* to use an embedded cubature based on the [`GrundmannMoeller`](@ref) rule of high order 13 and low order 11, you can do

```@example embedded-cubature
using HAdaptiveIntegration
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature

t = Triangle((0, 0), (1, 0), (0, 1))
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
ec = embedded_cubature(GrundmannMoeller{2}(13, 11))
I, E = integrate(f, t; rule = ec)
```

Which cubature rule is best depends on the function being integrated, as well as on the desired accuracy; as a rule of thumb, higher-order cubature rules will perform better for globally smooth functions `f` or higher accuracy requirements. Here is a short study on the number of function evaluations required to achieve a given accuracy for the default Radon-Laurie cubature and the `GrundmannMoeller` cubature rule above:

```@example embedded-cubature
const cc = Ref(0)  # a counter for the number of function evaluations
g = x -> (cc[] += 1; f(x))
rtol = 1e-2
cc[] = 0; integrate(g, t; rtol); counter_default = cc[]
cc[] = 0; integrate(g, t; rule = ec, rtol); counter_custom = cc[]
counter_default, counter_custom
```

For `rtol = 1e-2`, we see that the default cubature rule requires fewer function evaluations. However, decreasing `rtol` changes the balance:

```@example embedded-cubature
rtol = 1e-8
cc[] = 0; integrate(g, t; rtol); counter_default = cc[]
cc[] = 0; integrate(g, t; rule = ec, rtol); counter_custom = cc[]
counter_default, counter_custom
```

This example illustrates that testing is necessary to determine which cubature rule is best for your specific application!

!!! tip "Available embedded cubature formulas"
    The list of available embedded cubature formulas is:

    ```@example rule-name
    using HAdaptiveIntegration # hide
    not_rules = Set([ # hide
        "AbstractRule", # hide
        "Rule", # hide
        "EmbeddedCubature", # hide
        "TabulatedEmbeddedCubature", # hide
        "embedded_cubature", # hide
        "orders", # hide
    ]) # hide
    for name in map(String, names(HAdaptiveIntegration.Rule)) # hide
        if name ∉ not_rules # hide
            println(name) # hide
        end # hide
    end # hide
    ```

To add a custom embedded cubature for a given domain, you must write a constructor, *e.g.*, `my_custom_cubature(args...)` that returns a valid [`EmbeddedCubature`](@ref). See the file at [`Rule/triangle.jl`](https://github.com/zmoitier/HAdaptiveIntegration.jl/blob/main/src/Rule/triangle.jl) for examples. PRs with new schemes are more than welcome!

## Define custom subdivision strategies

The package uses a default subdivision strategy for the given `domain` by calling [`default_subdivision`](@ref). For example, by default triangles are subdivided into 4 smaller triangles by connecting the midpoints of the edges:

```@example default-subdivision
using HAdaptiveIntegration

t = Triangle((0, 0), (1, 0), (0, 1))
subdiv_algo = HAdaptiveIntegration.default_subdivision(t)
```

Here are the subdivided triangles:

```@example default-subdivision
subdiv_algo(t)
```

But it is also possible (and maybe desirable) to split the triangle into 2 smaller triangles instead. The following function accomplishes this:

```@example default-subdivision
using StaticArrays

function subdivide_triangle2(t::Triangle{T}) where {T}
    a, b, c = t.vertices
    bc = (b + c) / 2
    return (Triangle{T}(SVector(bc, a, b)), Triangle{T}(SVector(bc, c, a)))
end
subdivide_triangle2(t)  # hide
```

!!! warning
    Non-default subdivision strategies may affect convergence. The default (4-way split) typically requires fewer subdivisions than the 2-way split shown above.

Passing `subdivide_triangle2` as the `subdiv_algo` to `integrate` will use this instead of the default:

```@example default-subdivision
f = x -> 1 / (x[1]^2 + x[2]^2 + 1e-2)
I, E = integrate(f, t; subdiv_algo = subdivide_triangle2)
```

Which subdivision strategy is best depends on the function being integrated; for the example presented above, it turns out the default strategy is more efficient!
