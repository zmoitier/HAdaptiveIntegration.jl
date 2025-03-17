function _type_float(containers...)
    T = reduce(
        promote_type, mapreduce(typeof, promote_type, container) for container in containers
    )
    return float(T)
end

"""
    abstract type AbstractDomain{D,T<:Real}

Abstract type for domain's integration in `D` dimensions.

## Mandatory methods:
- [`map_from_reference`](@ref)
- [`abs_det_jac`](@ref)

## Useful (but non-mandatory) methods:
- [`reference_domain`](@ref)
- [`map_to_reference`](@ref)
"""
abstract type AbstractDomain{D,T<:Real} end

"""
    struct Orthotope{D,T} <: Domain{D,T}

Axes-aligned Orthotope in `D` dimensions, with value type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `low_corner::SVector{D,T}`: the low corner.
- `high_corner::SVector{D,T}`: the high corner.
"""
struct Orthotope{D,T} <: AbstractDomain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)
    orthotope(T::DataType, low_corner, high_corner)

Return an axes-aligned orthotope in `D` dimensions, with value type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function orthotope(T::DataType, low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end
function orthotope(low_corner, high_corner)
    return orthotope(_type_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    reference_orthotope(D::Int)
    reference_orthotope(T::DataType, D::Int)

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with value type `T`.
"""
function reference_orthotope(T::DataType, D::Int)
    return Orthotope(zeros(SVector{D,T}), ones(SVector{D,T}))
end
reference_orthotope(D::Int) = reference_orthotope(float(Int), D)

"""
    Segment{T} = Orthotope{1,T}
"""
const Segment{T} = Orthotope{1,T}

"""
    segment(xmin, xmax)
    segment(T::DataType, xmin, xmax)

Return a segment in 1 dimensions with value type `T` representing the interval
`[xmin, xmax]` with `xmin ≤ xmax`.
"""
function segment(T::DataType, xmin::R, xmax::S) where {R<:Real,S<:Real}
    return orthotope(T, (xmin,), (xmax,))
end
function segment(xmin::R, xmax::S) where {R<:Real,S<:Real}
    T = promote_type(typeof(xmin), typeof(xmax))
    return segment(float(T), xmin, xmax)
end

"""
    Rectangle{T} = Orthotope{2,T}
"""
const Rectangle{T} = Orthotope{2,T}

"""
    rectangle(low_corner, high_corner)
    rectangle(T::DataType, low_corner, high_corner)

Return an axes-aligned rectangle, with value type `T`, given by two 2d-points `low_corner`
and `high_corner` with. Note that, we must have `low_corner .≤ high_corner`.
"""
function rectangle(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner` must be 2d-vector."
    return orthotope(T, low_corner, high_corner)
end
function rectangle(low_corner, high_corner)
    return rectangle(_type_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    Cuboid{T} = Orthotope{3,T}
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)
    cuboid(T::DataType, low_corner, high_corner)

Return an axes-aligned cuboid, with value type `T`, given by two 3d-points `low_corner` and
`high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function cuboid(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner` must be 3d-vector."
    return orthotope(T, low_corner, high_corner)
end
function cuboid(low_corner, high_corner)
    return cuboid(_type_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    struct Simplex{D,T,N} <: Domain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of value type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.
"""
struct Simplex{D,T,N} <: AbstractDomain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{D,T,N}(vertices::SVector{D,T}...) where {D,T,N}
    return Simplex(SVector{N}(vertices...))
end

"""
    simplex(vertices...)
    simplex(T::DataType, vertices...)

Return a `D`-simplex with value type `T` from a collection of vertices. Note that all points
must have the same length `D` and there must be `N=D+1` points.
"""
function simplex(T::DataType, vertices...)
    N = length(vertices)
    D = N - 1
    @assert all(p -> length(p) == D, vertices)

    return Simplex(SVector{N}(SVector{D,T}.(vertices)))
end
simplex(vertices...) = simplex(_type_float(vertices...), vertices...)

"""
    reference_simplex(D::Int)
    reference_simplex(T::DataType, D::Int)

Return the reference `D`-dimensional simplex with value type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(T::DataType, D::Int)
    points = [zeros(SVector{D,T})]
    append!(points, setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)
    return Simplex(SVector{D + 1}(points))
end
reference_simplex(D::Int) = reference_simplex(float(Int), D)

"""
    Triangle{T} = Simplex{2,T,3}
"""
const Triangle{T} = Simplex{2,T,3}

"""
    triangle(a, b, c)
    triangle(T::DataType, a, b, c)

Return a triangle in 2 dimensions, with value type `T`, given by three 2d-points `a`, `b`,
and `c`.
"""
function triangle(T::DataType, a, b, c)
    @assert length(a) == length(b) == length(c) == 2 "`a`, `b`, and `c` must be 2d-vector."
    return simplex(T, a, b, c)
end
triangle(a, b, c) = triangle(_type_float(a, b, c), a, b, c)

"""
    Tetrahedron{T} = Simplex{3,T,4}
"""
const Tetrahedron{T} = Simplex{3,T,4}

"""
    tetrahedron(a, b, c, d)
    tetrahedron(T::DataType, a, b, c, d)

Return a tetrahedron in 3 dimensions, with value type `T`, given by four 3d-points `a`, `b`,
`c`, and `d`.
"""
function tetrahedron(T::DataType, a, b, c, d)
    @assert length(a) == length(b) == length(c) == length(d) == 3 "`a`, `b`, `c`, and `d`
    must be 3d-vector."
    return simplex(T, a, b, c, d)
end
tetrahedron(a, b, c, d) = tetrahedron(_type_float(a, b, c, d), a, b, c, d)

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
