"""
    struct Orthotope{D,T}

Axes-aligned Orthotope in `D` dimensions given by two points `low_corner` and `high_corner`.
"""
struct Orthotope{D,T} <: Domain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)

An axes-aligned orthotope in `D` dimensions given by two vectors `low_corner` and `high_corner`.
"""
function orthotope(low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a < b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    return Orthotope(SVector{D}(low_corner), SVector{D}(high_corner))
end

"""
    reference_orthotope(D, T=Float64)

Return the reference orthotope in `D` dimensions, representing `[0, 1]ᴰ`.
"""
function reference_orthotope(D::Int, T::DataType=Float64)
    @assert D ≥ 1 "D = $D must be greater than 1."

    low_corner = SVector{D,T}(zeros(T, D))
    high_corner = SVector{D,T}(ones(T, D))

    return Orthotope(low_corner, high_corner)
end

function map_from_reference(h::Orthotope{D,T}) where {D,T<:Real}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

function map_to_reference(h::Orthotope{D,T}) where {D,T<:Real}
    return u -> (u - h.low_corner) ./ (h.high_corner - h.low_corner)
end

function abs_det_jacobian(s::Orthotope{D,T}) where {D,T<:Real}
    return prod(s.high_corner - s.low_corner)
end

"""
A segment in 1 dimensions given by two value `xmin` and `xmax`.
"""
const Segment{T} = Orthotope{1,T}

function Segment(low_corner::SVector{1,T}, high_corner::SVector{1,T}) where {T<:Real}
    return Orthotope(low_corner, high_corner)
end

"""
    segment(xmin, xmax)

A segment in 1 dimensions representing `[xmin, xmax]`.
"""
function segment(xmin::T, xmax::T) where {T<:Real}
    return orthotope([xmin], [xmax])
end

"""
An axes-aligned rectangle given by two 2d-points `low_corner` and `high_corner`.
"""
const Rectangle{T} = Orthotope{2,T}

function Rectangle(low_corner::SVector{2,T}, high_corner::SVector{2,T}) where {T<:Real}
    return Orthotope(low_corner, high_corner)
end

"""
    rectangle(low_corner, high_corner)

An axes-aligned rectangle given by two 2d-points `low_corner` and `high_corner`.
"""
function rectangle(low_corner, high_corner)
    return orthotope(low_corner, high_corner)
end

"""
A axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
"""
const Cuboid{T} = Orthotope{3,T}

function Cuboid(low_corner::SVector{3,T}, high_corner::SVector{3,T}) where {T<:Real}
    return Orthotope(low_corner, high_corner)
end

"""
    cuboid(low_corner, high_corner)

A axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
"""
function cuboid(low_corner, high_corner)
    return orthotope(low_corner, high_corner)
end
