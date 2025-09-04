"""
    struct Segment{T} <: AbstractDomain{1,T}

Segment in 1 dimension, with element type `T`, given by `xmin` and `xmax`. Note that, we
must have `xmin ≤ xmax`.

## Fields:
- `xmin::T`: the low corner.
- `xmax::T`: the high corner.

## Invariants (**not** check at construction):
- `xmin ≤ xmax`

## Constructors:
- `Segment(xmin, xmax)`
- `Segment{T}(low_corner, high_corner)`
- `Segment(corners::SVector{2,SVector{1,T}})`
"""
struct Segment{T} <: AbstractDomain{1,T}
    xmin::T
    xmax::T

    function Segment{T}(xmin::T, xmax::T) where {T}
        return new{T}(xmin, xmax)
    end
end

function Segment{T}(xmin, xmax) where {T}
    @assert xmin < xmax "degenerate 1-dimensional Segment: must have `xmax > xmin`."
    return Segment{T}(T(xmin), T(xmax))
end

function Segment(xmin, xmax)
    @assert xmin < xmax "degenerate 1-dimensional Segment: must have `xmax > xmin`."
    T = float(promote_type(typeof(xmin), typeof(xmax)))
    return Segment{T}(T(xmin), T(xmax))
end

"""
    reference_segment(T=float(Int))

Return the reference 1-dimensional segment `[0, 1]` with element type `T`.
"""
function reference_segment((::Type{T})=float(Int)) where {T}
    return Segment(zero(T), one(T))
end

function map_from_reference(s::Segment{T}) where {T}
    return u -> s.xmin + u * (s.xmax - s.xmin)
end

function abs_det_jac(s::Segment{T}) where {T}
    return s.xmax - s.xmin
end

function map_to_reference(s::Segment{T}) where {T}
    diff = s.xmax - s.xmin
    @assert diff > √eps(float(T)) "degenerate 1-dimensional Segment: must have `high_corner .> low_corner`."

    return p -> (p - s.xmin) / diff
end

"""
    subdivide_segment(s::Segment)

Subdivide the 1-dimensional segment `s` into 2 smaller segments by splitting it at
its midpoint.
"""
function subdivide_segment(s::Segment{T}) where {T}
    m = (s.xmin + s.xmax) / 2
    return (Segment{T}(s.xmin, m), Segment{T}(m, s.xmax))
end
