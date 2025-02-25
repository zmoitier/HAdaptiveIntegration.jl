"""
    struct Orthotope{D,T} <: Domain{D,T}

Axes-aligned Orthotope in `D` dimensions given by two vectors `low_corner` and `high_corner`.
"""
struct Orthotope{D,T} <: Domain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner::T, high_corner::T) where {T}

An axes-aligned orthotope in `D` dimensions given by two vectors `low_corner` and `high_corner`.
"""
function orthotope(low_corner::T, high_corner::T) where {T}
    @assert (length(low_corner) == length(high_corner))
    @assert all(a < b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    return Orthotope(SVector{D}(low_corner), SVector{D}(high_corner))
end

"""
    reference_orthotope(D::Int, T::DataType)

Return the reference orthotope in `D` dimensions, representing `[0, 1]ᴰ`.
"""
function reference_orthotope(D::Int, T::DataType=Float64)
    @assert D ≥ 1 "D = $D must be greater than 1."

    low_corner = SVector{D,T}(zeros(T, D))
    high_corner = SVector{D,T}(ones(T, D))

    return Orthotope(low_corner, high_corner)
end

"""
    map_from_reference(h::Orthotope)::Function

Return an anonymous function that maps the reference orthotope to the physical orthotope `h`.
"""
function map_from_reference(h::Orthotope)::Function
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

"""
    map_to_reference(h::Orthotope)::Function

Return an anonymous function that maps the physical orthotope `h` to the reference orthotope.
"""
function map_to_reference(h::Orthotope)::Function
    return u -> (u - h.low_corner) ./ (h.high_corner - h.low_corner)
end

"""
    abs_jacobian_determinant(s::Orthotope)

The absolute value of the Jacobian's determinant of the map from the reference orthotope to the physical orthotope.
"""
function abs_jacobian_determinant(s::Orthotope)
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
    segment(xmin::T, xmax::T) where {T<:Real}

A segment in 1 dimensions representing `[xmin, xmax]`.
"""
function segment(xmin::T, xmax::T) where {T<:Real}
    return orthotope([xmin], [xmax])
end

"""
    reference_segment(D::Int, T::DataType)

Return the reference segment in 1 dimensions, representing `[0, 1]`.
"""
function reference_segment(T::DataType=Float64)
    return reference_orthotope(1, T)
end

"""
A rectangle in 2 dimensions given by two vectors `low_corner` and `high_corner`.
"""
const Rectangle{T} = Orthotope{2,T}

function Rectangle(low_corner::SVector{2,T}, high_corner::SVector{2,T}) where {T<:Real}
    return Orthotope(low_corner, high_corner)
end

"""
    rectangle(xlim::Tuple{T,T}, ylim::Tuple{T,T}) where {T<:Real}

A rectangle in 2 dimensions representing `xlim × ylim`.
"""
function rectangle(xlim::Tuple{T,T}, ylim::Tuple{T,T}) where {T<:Real}
    xmin, xmax = xlim
    ymin, ymax = ylim
    return orthotope(SVector{2,T}([xmin, ymin]), SVector{2,T}([xmax, ymax]))
end

"""
    reference_rectangle(D::Int, T::DataType)

Return the reference rectangle (square) in 2 dimensions, representing `[0, 1]²`.
"""
function reference_rectangle(T::DataType=Float64)
    return reference_orthotope(2, T)
end

"""
A cuboid in 3 dimensions given by two vectors `low_corner` and `high_corner`.
"""
const Cuboid{T} = Orthotope{3,T}

function Cuboid(low_corner::SVector{3,T}, high_corner::SVector{3,T}) where {T<:Real}
    return Orthotope(low_corner, high_corner)
end

"""
    cuboid(xlim::Tuple{T,T}, ylim::Tuple{T,T}, zlim::Tuple{T,T}) where {T<:Real}

A cuboid in 3 dimensions representing `xlim × ylim × zlim`.
"""
function cuboid(xlim::Tuple{T,T}, ylim::Tuple{T,T}, zlim::Tuple{T,T}) where {T<:Real}
    xmin, xmax = xlim
    ymin, ymax = ylim
    zmin, zmax = zlim
    return orthotope(SVector{3,T}([xmin, ymin, zmin]), SVector{3,T}([xmax, ymax, zmax]))
end

"""
    reference_cuboid(D::Int, T::DataType)

Return the reference cuboid (cube) in 3 dimensions, representing `[0, 1]³`.
"""
function reference_cuboid(T::DataType=Float64)
    return reference_orthotope(3, T)
end
