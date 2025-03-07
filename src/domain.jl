"""
    abstract type Domain{D,T<:Real}

Abstract type for integration domains in `D` dimensions with type `T`.
"""
abstract type Domain{D,T<:Real} end

"""
    dimension(domain::Domain{D,T}) where {D,T<:Real}

Return the dimension `D` of `domain`.
"""
function dimension(::Domain{D,T}) where {D,T<:Real}
    return D
end

"""
    value_type(domain::Domain{D,T}) where {D,T<:Real}

Return the type `T` of `domain`.
"""
function value_type(::Domain{D,T}) where {D,T<:Real}
    return T
end

"""
    struct Orthotope{D,T} <: Domain{D,T}

Axes-aligned Orthotope in `D` dimensions given by two points `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.
"""
struct Orthotope{D,T} <: Domain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)

An axes-aligned orthotope in `D` dimensions given by two vectors `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.
"""
function orthotope(low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    T = promote_type(eltype(low_corner), eltype(high_corner))
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end

"""
    Segment{T} = Orthotope{1,T}
"""
const Segment{T} = Orthotope{1,T}

"""
    segment(xmin::T, xmax::T) where {T<:Real}

A segment in 1 dimensions representing the interval `[xmin, xmax]` with `xmin ≤ xmax`.
"""
function segment(xmin::T, xmax::T) where {T<:Real}
    return orthotope(SVector(xmin), SVector(xmax))
end

"""
    Rectangle{T} = Orthotope{2,T}
"""
const Rectangle{T} = Orthotope{2,T}

"""
    rectangle(low_corner, high_corner)

An axes-aligned rectangle given by two 2d-points `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.
"""
function rectangle(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner` must be 2d-vector."
    return orthotope(SVector{2}(low_corner), SVector{2}(high_corner))
end

"""
    Cuboid{T} = Orthotope{3,T}
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)

A axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.
"""
function cuboid(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner` must be 3d-vector."
    return orthotope(SVector{3}(low_corner), SVector{3}(high_corner))
end

"""
    struct Simplex{D,T,N} <: Domain{D,T}

A simplex in `D` dimensions with `N=D+1` points of element type `T`.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    points::SVector{N,SVector{D,T}}
end

function Simplex{D,T,N}(points::SVector{D,T}...) where {D,T<:Real,N}
    return Simplex(SVector{N}(points...))
end

"""
    simplex(points...)

A simplex in `D` dimensions with `N=D+1` points of element type `T`.
"""
function simplex(points...)
    N = length(points)
    D = N - 1
    @assert all(p -> length(p) == D, points)

    return Simplex(SVector{N}(SVector{D}.(points)))
end

"""
    Triangle{T} = Simplex{2,T,3}
"""
const Triangle{T} = Simplex{2,T,3}

"""
    triangle(a, b, c)

A triangle in 2 dimensions given by three 2d-points `a`, `b`, and `c`.
"""
function triangle(a, b, c)
    @assert length(a) == length(b) == length(c) == 2 "`a`, `b`, and `c` must be 2d-vector."
    return simplex(a, b, c)
end

"""
    Tetrahedron{T} = Simplex{3,T,4}
"""
const Tetrahedron{T} = Simplex{3,T,4}

"""
    tetrahedron(a, b, c, d)

A tetrahedron in 3 dimensions given by four 3d-points `a`, `b`, `c`, and `d`.
"""
function tetrahedron(a, b, c, d)
    @assert length(a) == length(b) == length(c) == length(d) == 2 "`a`, `b`, `c`, and `d` must be 2d-vector."
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
