"""
    abstract type AbstractDomain{D,T}

Abstract type for integration domains' in `D` dimensions with element type `T`.

## Type Parameters:
- `D`: dimension of the domain.
- `T`: element type of the domain.

## Mandatory methods:
- [`map_from_reference`](@ref)
- [`abs_det_jac`](@ref)

## Useful (but non-mandatory) methods:
- [`reference_domain`](@ref)
- [`map_to_reference`](@ref)
"""
abstract type AbstractDomain{D,T} end

"""
    dimension(::Type{DOM}) where {DOM<:AbstractDomain{D}}

Return the dimension `D` of the given domain `DOM`.
"""
dimension(::Type{AbstractDomain{D,T}}) where {D,T} = D
dimension(::Type{AbstractDomain{D}}) where {D} = D
dimension(::Type{DOM}) where {DOM<:AbstractDomain} = dimension(supertype(DOM))

"""
    struct Orthotope{D,T} <: AbstractDomain{D,T}

Axes-aligned Orthotope in `D` dimensions, with element type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `low_corner::SVector{D,T}`: the low corner.
- `high_corner::SVector{D,T}`: the high corner.

## Invariants (**not** check at construction):
- `low_corner .≤ high_corner`
"""
struct Orthotope{D,T} <: AbstractDomain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)
    orthotope(T::DataType, low_corner, high_corner)

Return an axes-aligned orthotope in `D` dimensions, with element type `T`, given by two
points `low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function orthotope(T::DataType, low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end
function orthotope(low_corner, high_corner)
    return orthotope(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    reference_orthotope(D::Int)
    reference_orthotope(T::DataType, D::Int)

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with element type `T`.
"""
function reference_orthotope(T::DataType, D::Int)
    return Orthotope(zeros(SVector{D,T}), ones(SVector{D,T}))
end
reference_orthotope(D::Int) = reference_orthotope(float(Int), D)

"""
    Segment{T} = Orthotope{1,T}

Alias for a 1-dimensional segment of element type `T`.
"""
const Segment{T} = Orthotope{1,T}

"""
    segment(xmin, xmax)
    segment(T::DataType, xmin, xmax)

Return a segment in 1 dimensions with element type `T` representing the interval
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

Alias for a 2-dimensional rectangle of element type `T`.
"""
const Rectangle{T} = Orthotope{2,T}

"""
    rectangle(low_corner, high_corner)
    rectangle(T::DataType, low_corner, high_corner)

Return an axes-aligned rectangle, with element type `T`, given by two 2d-points `low_corner`
and `high_corner` with. Note that, we must have `low_corner .≤ high_corner`.
"""
function rectangle(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner` must be 2d-vector."
    return orthotope(T, low_corner, high_corner)
end
function rectangle(low_corner, high_corner)
    return rectangle(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    Cuboid{T} = Orthotope{3,T}

Alias for a 3-dimensional cuboid of value type `T`.
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)
    cuboid(T::DataType, low_corner, high_corner)

Return an axes-aligned cuboid, with element type `T`, given by two 3d-points `low_corner`
and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function cuboid(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner` must be 3d-vector."
    return orthotope(T, low_corner, high_corner)
end
function cuboid(low_corner, high_corner)
    return cuboid(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    struct Simplex{D,N,T} <: AbstractDomain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of element type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.

## Invariants (**not** check at construction):
- `N = D+1`
"""
struct Simplex{D,N,T} <: AbstractDomain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{D,N,T}(vertices::Vararg{SVector{D,T},N}) where {D,N,T}
    return Simplex(SVector{N}(vertices...))
end

"""
    simplex(vertices...)
    simplex(T::DataType, vertices...)

Return a `D`-simplex with element type `T` from a collection of `N=D+1` vertices. Note that
all vertices must have the same length `D`.
"""
function simplex(T::DataType, vertices...)
    @assert allequal(length, vertices) "all `vertices` must have the same length."
    D = length(first(vertices))

    N = D + 1
    @assert length(vertices) == N "Expected $N vertices, but got $(length(vertices))."

    return Simplex(SVector{N}(SVector{D,T}.(vertices)))
end
simplex(vertices...) = simplex(promote_to_float(vertices...), vertices...)

"""
    reference_simplex(D::Int)
    reference_simplex(T::DataType, D::Int)

Return the reference `D`-dimensional simplex with element type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(T::DataType, D::Int)
    vertices = [
        zeros(SVector{D,T}), collect(setindex(zeros(SVector{D,T}), 1, i) for i in 1:D)...
    ]
    return Simplex(SVector{D + 1}(vertices))
end
reference_simplex(D::Int) = reference_simplex(float(Int), D)

"""
    Triangle{T} = Simplex{2,3,T}

Alias for a 2-dimensional triangle with 3 vertices of value type `T`.
"""
const Triangle{T} = Simplex{2,3,T}

"""
    triangle(a, b, c)
    triangle(T::DataType, a, b, c)

Return a triangle in 2 dimensions, with element type `T`, given by three 2d-points `a`, `b`,
and `c`.
"""
function triangle(T::DataType, a, b, c)
    @assert length(a) == length(b) == length(c) == 2 "`a`, `b`, and `c` must be 2d-vector."
    return simplex(T, a, b, c)
end
triangle(a, b, c) = triangle(promote_to_float(a, b, c), a, b, c)

"""
    Tetrahedron{T} = Simplex{3,4,T}

Alias for a 3-dimensional tetrahedron with 4 vertices of element type `T`.
"""
const Tetrahedron{T} = Simplex{3,4,T}

"""
    tetrahedron(a, b, c, d)
    tetrahedron(T::DataType, a, b, c, d)

Return a tetrahedron in 3 dimensions, with element type `T`, given by four 3d-points `a`,
`b`, `c`, and `d`.
"""
function tetrahedron(T::DataType, a, b, c, d)
    @assert length(a) == length(b) == length(c) == length(d) == 3 "`a`, `b`, `c`, and `d` must be 3d-vector."
    return simplex(T, a, b, c, d)
end
tetrahedron(a, b, c, d) = tetrahedron(promote_to_float(a, b, c, d), a, b, c, d)

"""
    map_from_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the reference domain to the physical domain `domain`.
"""
function map_from_reference(h::Orthotope{D,T}) where {D,T}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

function map_from_reference(s::Simplex{D,N,T}) where {D,N,T}
    return u -> begin
        v = (1 - sum(u)) * s.vertices[1]
        for i in 2:N
            v += u[i - 1] * s.vertices[i]
        end
        return v
    end
end

"""
    abs_det_jac(domain::DOM) where {DOM<:AbstractDomain}

Return the absolute value of the Jacobian's determinant of the map from the reference domain
to the physical domain `domain`.
"""
function abs_det_jac(h::Orthotope{D,T}) where {D,T}
    return prod(h.high_corner - h.low_corner)
end

function abs_det_jac(s::Simplex{D,N,T}) where {D,N,T}
    vertices = s.vertices
    jacobian_matrix = hcat(ntuple(i -> vertices[i + 1] - vertices[1], D)...)
    return abs(det(jacobian_matrix))
end

"""
    map_to_reference(domain::DOM) where {DOM<:AbstractDomain}

Return an anonymous function that maps the physical domain `domain` to the reference domain.

## Constraints:
- For `Orthotope`, must have `high_corner .> low_corner`.
- For `Simplex{D}`, the vertices must form a valid `D`-dimensional simplex with non-zero
  volume).
"""
function map_to_reference(h::Orthotope{D,T}) where {D,T}
    diff = h.high_corner - h.low_corner
    @assert all(x -> x > √eps(float(T)), diff) "degenerate $D-dimensional Orthotope: must have `high_corner .> low_corner`."

    return p -> (p - h.low_corner) ./ diff
end

function map_to_reference(s::Simplex{D,N,T}) where {D,N,T}
    v = s.vertices
    jacobian_matrix = hcat(ntuple(i -> v[i + 1] - v[1], D)...)
    @assert !isapprox(det(jacobian_matrix), 0; atol=√eps(float(T))) "degenerate $D-dimensional Simplex: the Jacobian matrix is not invertible."

    M = inv(jacobian_matrix)
    return u -> M * (u - v[1])
end

"""
    reference_domain(::Type{<:AbstractDomain})

Return the reference domain for the given domain type.
"""
reference_domain(::Type{Orthotope{D,T}}) where {D,T} = reference_orthotope(T, D)
reference_domain(::Type{Orthotope{D}}) where {D} = reference_orthotope(float(Int), D)

reference_domain(::Type{Simplex{D,N,T}}) where {D,N,T} = reference_simplex(T, D)
reference_domain(::Type{Simplex{D,N}}) where {D,N} = reference_simplex(float(Int), D)
reference_domain(::Type{Simplex{D}}) where {D} = reference_simplex(float(Int), D)
