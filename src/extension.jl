"""
    increase_precision(
        tec::TabulatedEmbeddedCubature,
        ::Type{T};
        x_atol=10 * eps(T),
        f_atol=10 * eps(T),
        maxiter::Int=16,
    )

Increase the precision of a `TabulatedEmbeddedCubature` to match the precision of the
floating-point type `T`. This method is provided by an extension, so `ForwardDiff` must be
loaded.

## Arguments
- `tec::TabulatedEmbeddedCubature`: the tabulated embedded cubature to increase the 
  precision of.
- `::Type{T}`: target floating-point type.

## Optional arguments
- `x_atol=10 * eps(T)`: the absolute tolerance for the change in the variables (nodes and
  weights) between Newton iterations.
- `f_atol=10 * eps(T)`: the absolute tolerance for the change in the function values
  (integral constraints) between Newton iterations.
- `maxiter::Int=16`: the maximum number of Newton iterations to perform.
"""
function increase_precision end
