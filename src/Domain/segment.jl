"""
    struct Segment{T} <: AbstractDomain{1,T}

A segment in 1D with element type `T`, given by two end points `a` and `b`.

## Fields:
- `a::T`
- `b::T`
"""
struct Segment{T} <: AbstractDomain{1,T}
    a::T
    b::T
end

"""
    reference_segment(T::DataType)

Return the reference `1`-dimensional segment `[0, 1]` with element type `T`.
"""
function reference_segment(T::DataType)
    return Segment{T}(zero(T), one(T))
end

function map_from_reference(s::Segment{T}) where {T}
    return u -> s.a .+ u .* (s.b - s.a)
end

function abs_det_jac(s::Segment{T}) where {T}
    return abs(s.b - s.a)
end

function map_to_reference(s::Segment{T}) where {T}
    diff = s.b - s.a
    @assert abs(diff) > √eps(float(T)) "degenerate Segment: `a ≈ b`."

    return p -> (p .- s.a) ./ diff
end

"""
    subdivide_segment(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment(s::Segment{T}) where {T}
    a, b = s.a, s.b
    m = (a + b) / 2
    return (Segment{T}(a, m), Segment{T}(m, b))
end
