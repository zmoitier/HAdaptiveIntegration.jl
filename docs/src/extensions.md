```@meta
CurrentModule = HAdaptiveIntegration
```

# [Extensions](@id extensions)

## `IncreasePrecisionExt`

This extension allows increasing the precision of tabulated embedded cubature. To use this
extension you must add the package
[`ForwardDiff`](https://juliadiff.org/ForwardDiff.jl/stable/).

```@example increase_precision
using HAdaptiveIntegration
using ForwardDiff
```

For example if we want to increase the precision of the rule
[`HAdaptiveIntegration.SQUARE_CH21`](@ref), we do:

```@example increase_precision
tec0 = HAdaptiveIntegration.SQUARE_CH21
tec0.nodes
```

```@example increase_precision
tec1 = HAdaptiveIntegration.increase_precision(
    tec0, BigFloat; x_atol=big"1e-64", f_atol=big"1e-64"
)
display(tec1) # hide
```

```@example increase_precision
tec1.nodes
```

## Complete workflow: Arbitrary precision integration

Now let's use the increased precision rule in an actual integration. First define a function
to integrate over a domain:

```@example increase_precision
domain = Rectangle((big"0", big"0"), (big"1", big"1"))
f = x -> x[1]^2 + x[2]^2
nothing # hide
```

Create an embedded cubature with the high-precision rule and a lower-order pair:

```@example increase_precision
using HAdaptiveIntegration.Rule: embedded_cubature
ec_big = embedded_cubature(tec1, BigFloat)
nothing # hide
```

Now integrate with arbitrary precision:

```@example increase_precision
I, E = HAdaptiveIntegration.integrate(f, domain; rule = ec_big, rtol = big"1e-64")
println("I = $I") # hide
println("E = $E") # hide
```

The result is now computed with arbitrary precision (BigFloat), enabling high-precision
numerical integration when needed for sensitive applications or validation studies.
