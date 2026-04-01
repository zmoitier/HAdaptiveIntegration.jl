"""
    struct Segment{T} <: AbstractDomain{1,T}

One-dimensional segment with element type `T`, defined by `xmin` and `xmax`.

## Fields:
- `xmin::T`: lower endpoint.
- `xmax::T`: upper endpoint.

## Invariants (checked by outer constructors):
- `xmin < xmax`

## Constructors:
- `Segment(xmin, xmax)`
- `Segment{T}(xmin, xmax)`
"""
struct Segment{T} <: AbstractDomain{1, T}
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
function reference_segment((::Type{T}) = float(Int)) where {T}
    return Segment{T}(zero(T), oneunit(T))
end

function map_from_reference(s::Segment{T}) where {T}
    diff = s.xmax - s.xmin
    return (u -> s.xmin .+ u .* diff, diff)
end

function map_to_reference(s::Segment{T}) where {T}
    diff = s.xmax - s.xmin
    @assert diff > √eps(float(T)) "degenerate 1-dimensional Segment: must have \
`high_corner .> low_corner`."

    return p -> (p .- s.xmin) ./ diff
end

"""
    subdivide_segment(s::Segment)

Subdivide `s` into two smaller segments by splitting at its midpoint.
"""
function subdivide_segment(s::Segment{T}) where {T}
    m = (s.xmin + s.xmax) / 2
    return (Segment{T}(s.xmin, m), Segment{T}(m, s.xmax))
end
