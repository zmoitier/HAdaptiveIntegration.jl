"""
    increase_precision(
        tec::TabulatedEmbeddedCubature,
        ::Type{T};
        x_atol=10 * eps(T),
        f_atol=10 * eps(T),
        maxiter::Int=16,
    )

Increase the precision of a `TabulatedEmbeddedCubature` to match the precision of the
floating point type `T`. This is an extension, you need `using ForwardDiff` for using it.

## Arguments
- `tec::TabulatedEmbeddedCubature`: the tabulated embedded cubature to increase the
  precision of.
- `::Type{T}`: the floating point type to increase the precision to.

## Optional arguments
- `x_atol=10 * eps(T)`: the absolute tolerance for the change in the variables (nodes and
  weights) between Newton iterations.
- `f_atol=10 * eps(T)`: the absolute tolerance for the change in the function values
  (integral constraints) between Newton iterations.
- `maxiter::Int=16`: the maximum number of Newton iterations to perform.
"""
function increase_precision end
