"""
    abstract type Domain{D,T<:Real}

Abstract type for integration on domains in `D` dimensions with value type `T`.
"""
abstract type Domain{D,T<:Real} end

"""
    struct Orthotope{D,T} <: Domain{D,T}

Axes-aligned Orthotope in `D` dimensions given by two points `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `low_corner::SVector{D,T}`: the low corner.
- `high_corner::SVector{D,T}`: the high corner.
"""
struct Orthotope{D,T} <: Domain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)

Return an axes-aligned orthotope in `D` dimensions given by two points `low_corner` and
`high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function orthotope(low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    T = promote_type(eltype(low_corner), eltype(high_corner))
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end

"""
    reference_orthotope(D::Int, T::DataType=Float64)

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with value type `T`.
"""
function reference_orthotope(D::Int, T::DataType=Float64)
    return Orthotope(zeros(SVector{D,T}), ones(SVector{D,T}))
end

"""
    Segment{T} = Orthotope{1,T}
"""
const Segment{T} = Orthotope{1,T}

"""
    segment(xmin, xmax)

Return a segment in 1 dimensions representing the interval `[xmin, xmax]` with
`xmin ≤ xmax`.
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

Return an axes-aligned rectangle given by two 2d-points `low_corner` and `high_corner`. Note
that, we must have `low_corner .≤ high_corner`.
"""
function rectangle(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner`
    must be 2d-vector."
    return orthotope(SVector{2}(low_corner), SVector{2}(high_corner))
end

"""
    Cuboid{T} = Orthotope{3,T}
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)

Return an axes-aligned cuboid given by two 3d-points `low_corner` and `high_corner`.
Note that, we must have `low_corner .≤ high_corner`.
"""
function cuboid(low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner`
    must be 3d-vector."
    return orthotope(SVector{3}(low_corner), SVector{3}(high_corner))
end

"""
    struct Simplex{D,T,N} <: Domain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of value type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{D,T,N}(vertices::SVector{D,T}...) where {D,T,N}
    return Simplex(SVector{N}(vertices...))
end

"""
    simplex(vertices...)

Return a `D`-simplex from a collection of vertices.
Note that all points must have the same length `D` and there must be `N=D+1` points.
"""
function simplex(vertices...)
    N = length(vertices)
    D = N - 1
    @assert all(p -> length(p) == D, vertices)

    return Simplex(SVector{N}(SVector{D}.(vertices)))
end

"""
    reference_simplex(D::Int, T::DataType=Float64)

Return the reference `D`-dimensional simplex with value type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(D::Int, T::DataType=Float64)
    points = [zeros(SVector{D,T})]
    append!(points, setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)
    return Simplex(SVector{D + 1}(points))
end

"""
    Triangle{T} = Simplex{2,T,3}
"""
const Triangle{T} = Simplex{2,T,3}

"""
    triangle(a, b, c)

Return a triangle in 2 dimensions given by three 2d-points `a`, `b`, and `c`.
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

Return a tetrahedron in 3 dimensions given by four 3d-points `a`, `b`, `c`, and `d`.
"""
function tetrahedron(a, b, c, d)
    @assert length(a) == length(b) == length(c) == length(d) == 3 "`a`, `b`, `c`, and `d`
    must be 3d-vector."
    return simplex(a, b, c, d)
end

"""
    map_from_reference(domain::Domain)

Return an anonymous function that maps the reference domain to the physical domain `domain`.
"""
function map_from_reference(h::Orthotope{D,T}) where {D,T}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

function map_from_reference(s::Simplex{D,T,N}) where {D,T,N}
    return u -> begin
        v = (1 - sum(u)) * s.vertices[1]
        for i in 2:N
            v += u[i - 1] * s.vertices[i]
        end
        return v
    end
end

"""
    map_to_reference(domain::Domain)

Return an anonymous function that maps the physical domain `domain` to the reference domain.
"""
function map_to_reference(h::Orthotope{D,T}) where {D,T}
    return u -> (u - h.low_corner) ./ (h.high_corner - h.low_corner)
end

function map_to_reference(s::Simplex{D,T,N}) where {D,T,N}
    v = s.vertices
    M = inv(hcat(ntuple(i -> v[i + 1] - v[1], D)...))
    return u -> M * (u - v[1])
end

"""
    abs_det_jac(domain::Domain)

Return the absolute value of the Jacobian's determinant of the map from the reference domain
to the physical domain `domain`.
"""
function abs_det_jac(h::Orthotope{D,T}) where {D,T}
    return prod(h.high_corner - h.low_corner)
end

function abs_det_jac(s::Simplex{D,T,N}) where {D,T,N}
    v = s.vertices
    mat = hcat(ntuple(i -> v[i + 1] - v[1], D)...)
    return abs(det(mat))
end
