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
