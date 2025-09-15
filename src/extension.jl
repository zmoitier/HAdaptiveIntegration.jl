"""
    increase_precision(
        tec::TabulatedEmbeddedCubature, T, options=Optim.Options()
    )::EmbeddedCubature

Increase the precision of a tabulated embedded cubature rule `tec` to be suitable for type `T`.

## Arguments
- `tec::TabulatedEmbeddedCubature`: The tabulated embedded cubature.
- `T`: The target type for the increased precision.

## Optional arguments
- `options=Optim.Options()`: Optional optimization settings, the full parameters list can be
   found here:
   https://julianlsolvers.github.io/Optim.jl/stable/user/config/#General-Options
"""
function increase_precision end
