"""
    abstract type Domain{D,T<:Real}

Abstract type for integration domains in `D` dimensions with type `T<:Real`.
"""
abstract type Domain{D,T<:Real} end

dimension(::Domain{D,T}) where {D,T<:Real} = D
value_type(::Domain{D,T}) where {D,T<:Real} = T

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
A segment in 1 dimensions given by two value `xmin` and `xmax`.
"""
const Segment{T} = Orthotope{1,T}

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

"""
    rectangle(low_corner, high_corner)

An axes-aligned rectangle given by two 2d-points `low_corner` and `high_corner`.
"""
function rectangle(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "must be 2d-vector."
    return orthotope(low_corner, high_corner)
end

"""
A axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)

A axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
"""
function cuboid(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "must be 3d-vector."
    return orthotope(low_corner, high_corner)
end

"""
    struct Simplex{D,T,N}

A simplex in `D` dimensions with N=D+1 points of type `T`.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    points::SVector{N,SVector{D,T}}
end

function Simplex{D,T,N}(points::SVector{D,T}...) where {D,T<:Real,N}
    return Simplex(SVector{N}(points...))
end

"""
    simplex(points...)

A simplex in `D` dimensions with N=D+1 points of type `T`.
"""
function simplex(points...)
    N = length(points)
    D = N - 1
    @assert all(pt -> length(pt) == D, points)

    return Simplex(SVector{N}(SVector{D}.(points)))
end

"""
A triangle in 2 dimensions with 3 vertices of type `T`.
"""
const Triangle{T} = Simplex{2,T,3}

"""
    triangle(a, b, c)

A triangle in 2 dimensions given by the 2d-points `a`, `b`, and `c`.
"""
function triangle(a, b, c)
    return simplex(a, b, c)
end

"""
A tetrahedron in 3 dimensions with 4 points of type `T`.
"""
const Tetrahedron{T} = Simplex{3,T,4}

"""
    tetrahedron(a, b, c, d)

A tetrahedron in 3 dimensions given by the 3d-points `a`, `b`, `c`, and `d`.
"""
function tetrahedron(a, b, c, d)
    return simplex(a, b, c, d)
end

const LIST_DOMAIN = [
    "dimension 1" => ["segment"],
    "dimension 2" => ["rectangle", "triangle"],
    "dimension 3" => ["cuboid", "tetrahedron"],
    "dimension D" => ["orthotope", "simplex"],
]

"""
    reference_domain(domain_type::DataType)

Return the reference domain for `domain_type`.

Return the reference D-simplex given by the vertices `(0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1)`.

Return the reference orthotope in `D` dimensions, representing `[0, 1]ᴰ`.
"""
function reference_domain(domain_type::DataType)
    @assert (domain_type <: Domain) "only subtype of `Domain` are allowed."
    D = domain_type.parameters[1]
    T = domain_type.parameters[2]
    @assert D ≥ 1 "D = $D must be greater than 1."

    if domain_type <: Simplex
        N = D + 1
        points = [zeros(T, D) for _ in 1:N]
        for i in 1:D
            points[i + 1][i] = 1
        end
        return Simplex(SVector{N,SVector{D,T}}(SVector{D,T}.(points)))
    end

    if domain_type <: Orthotope
        low_corner = SVector{D,T}(zeros(T, D))
        high_corner = SVector{D,T}(ones(T, D))
        return Orthotope(low_corner, high_corner)
    end

    throw("no reference domain implemented for $domain_type.")
end

"""
    map_from_reference(d::Domain)::Function

Return an anonymous function that maps the reference domain to the physical domain `d`.
"""
function map_from_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is not implemented for type $(typeof(d))."
end

function map_from_reference(h::Orthotope{D,T}) where {D,T<:Real}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

function map_from_reference(s::Simplex{D,T,N}) where {D,T<:Real,N}
    # return u -> (1 - sum(u)) * s.points[1] + sum(u[i - 1] * s.points[i] for i in 2:N)
    # NOTE: for loop is faster than the above expression
    u -> begin
        v = (1 - sum(u)) * s.points[1]
        for i in 2:N
            v += u[i - 1] * s.points[i]
        end
        return v
    end
end

"""
    map_to_reference(d::Domain)::Function

Return an anonymous function that maps the physical domain `d` to the reference domain.
"""
function map_to_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_to_reference` is not implemented for type $(typeof(d))."
end

function map_to_reference(h::Orthotope{D,T}) where {D,T<:Real}
    return u -> (u - h.low_corner) ./ (h.high_corner - h.low_corner)
end

function map_to_reference(s::Simplex{D,T,N}) where {D,T<:Real,N}
    v = s.points
    M = inv(SMatrix{D,D}(reinterpret(reshape, T, [x - v[1] for x in v[2:N]])))
    return u -> M * (u - v[1])
end

"""
    abs_det_jacobian(d::Domain)

The absolute value of the Jacobian's determinant of the map from the reference domain to the physical domain `d`.
"""
function abs_det_jacobian(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is not implemented for $(typeof(d))."
end

function abs_det_jacobian(s::Orthotope{D,T}) where {D,T<:Real}
    return prod(s.high_corner - s.low_corner)
end

function abs_det_jacobian(s::Simplex{D,T,N}) where {D,T<:Real,N}
    v = s.points
    mat = hcat(ntuple(i -> v[i + 1] - v[1], D)...)
    return abs(det(mat))
end
