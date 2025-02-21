"""
    struct Orthotope{D,T} <: Domain{D,T}

Axes-aligned Orthotope in `D` dimensions given by two vectors `low_corner` and `high_corner`.
"""
struct Orthotope{D,T} <: Domain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

function Orthotope(lc::T, hc::T) where {T}
    @assert (length(lc) == length(hc))
    @assert all(a < b for (a, b) in zip(lc, hc))

    D = length(lc)
    return Orthotope(SVector{D}(lc), SVector{D}(hc))
end

"""
    reference_orthotope(D::Int, T::DataType)

Return the reference orthotope in `D` dimensions, given by `[0, 1]^D`.
"""
function reference_orthotope(D::Int, T::DataType=Float64)
    @assert D ≥ 1 "D = $D must be greater than 1."
    lc = @SVector zeros(T, D)
    hc = @SVector ones(T, D)
    return Orthotope(lc, hc)
end

"""
    map_from_reference(h::Orthotope)::Function

Return an anonymous function that maps the reference orthotope to the physical orthotope `h`.
"""
function map_from_reference(h::Orthotope)::Function
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

"""
    abs_jacobian_determinant(s::Orthotope)

The determinant of the Jacobian of the map from the reference orthotope to the physical orthotope.
"""
function abs_jacobian_determinant(s::Orthotope)
    return prod(s.high_corner - s.low_corner)
end

const Segment{T} = Orthotope{1,T}
const Rectangle{T} = Orthotope{2,T}
const Cuboid{T} = Orthotope{3,T}

# default types
Segment(args...) = Segment{Float64}(args...)
Rectangle(args...) = Rectangle{Float64}(args...)
Cuboid(args...) = Cuboid{Float64}(args...)

"""
    reference_segment(D::Int, T::DataType)

Return the reference orthotope in 2 dimensions, given by `[0, 1]^2`.
"""
function reference_segment(T::DataType=Float64)
    return reference_orthotope(1, T)
end

"""
    reference_rectangle(D::Int, T::DataType)

Return the reference rectangle (square) in 2 dimensions, given by `[0, 1]^2`.
"""
function reference_rectangle(T::DataType=Float64)
    return reference_orthotope(2, T)
end

"""
    reference_cuboid(D::Int, T::DataType)

Return the reference cuboid (cube) in 3 dimensions, given by `[0, 1]^3`.
"""
function reference_cuboid(T::DataType=Float64)
    return reference_orthotope(3, T)
end
