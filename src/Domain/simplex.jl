"""
    struct Simplex{D,T,N} <: AbstractDomain{D,T}

A simplex in `D` dimensions with `N=D+1` vertices of element type `T`.

## Fields:
- `vertices::SVector{N,SVector{D,T}}`: vertices of the simplex.

## Invariants (**not** check at construction):
- `N = D+1`

## Constructors:
- `Simplex(vertices...)`
- `Simplex{T}(vertices...)`
- `Simplex(vertices::SVector{N,SVector{D,T}})`
"""
struct Simplex{D,T,N} <: AbstractDomain{D,T}
    vertices::SVector{N,SVector{D,T}}
end

function Simplex{T}(vertices...) where {T}
    N = length(vertices)
    D = N - 1

    for vertex in vertices
        @assert length(vertex) == D "$vertex must have length $D."
    end

    return Simplex(SVector{N}(SVector{D,T}.(vertices)...))
end

function Simplex(vertices...)
    N = length(vertices)
    D = N - 1

    for vertex in vertices
        @assert length(vertex) == D "$vertex must have length $D."
    end

    return Simplex(SVector{N}(float.(SVector{D}.(vertices))...))
end

"""
    reference_simplex(D::Int, T=float(Int))

Return the reference `D`-dimensional simplex with element type `T`, which is the convex hull
of the `N=D+1` points `(0,...,0)`, `(1,0,...,0)`, `(0,1,0,...,0)`, ..., `(0,...,0,1)`.
"""
function reference_simplex(::Val{D}, (::Type{T})=float(Int)) where {D,T}
    vertices = ntuple(D + 1) do i
        i == 1 ? zeros(SVector{D,T}) : setindex(zeros(SVector{D,T}), 1, i - 1)
    end
    return Simplex(SVector{D + 1}(vertices))
end

function map_from_reference(s::Simplex{D,T,N}) where {D,T,N}
    vertices = s.vertices
    jacobian_matrix = hcat(ntuple(i -> vertices[i + 1] - vertices[1], D)...)
    return (u -> vertices[1] + jacobian_matrix * u, abs(det(jacobian_matrix)))
end

function map_to_reference(s::Simplex{D,T,N}) where {D,T,N}
    v = s.vertices
    jacobian_matrix = hcat(ntuple(i -> v[i + 1] - v[1], D)...)
    @assert !isapprox(det(jacobian_matrix), 0; atol=(√eps(float(T)))) "degenerate \
$D-dimensional Simplex: the Jacobian matrix is not invertible."

    M = inv(jacobian_matrix)
    return u -> M * (u - v[1])
end

"""
    subdivide_simplex_scheme(::Val{D}) where {D}

Return the color schemes for subdividing a `D`-simplex into 2ᴰ simplices by using the
`SimpleS` algorithm.
"""
@generated function subdivide_simplex_scheme(::Val{D}) where {D}
    N = D + 1

    color_schemes = ntuple(Val(2^D)) do k
        bits = k - 1
        cs = MMatrix{2,N,Int}(undef)

        cs[1, 1] = 1
        for j in 1:D
            cs[1, j + 1] = cs[1, j] + iseven(bits >> (j - 1))
        end

        cs[2, 1] = cs[1, N]
        for j in 1:D
            cs[2, j + 1] = cs[2, j] + isodd(bits >> (j - 1))
        end

        return SMatrix{2,N,Int}(cs)
    end

    return :($color_schemes)
end

"""
    subdivide_simplex(s::Simplex{D,T,N}) where {D,T,N}

Subdivide a `D`-simplex into 2ᴰ simplices by using the `SimpleS` algorithm.

Implements the `SimpleS` algorithm in [Algorithm 860: SimpleS -- an extension of
Freudenthal's simplex subdivision](https://doi.org/10.1145/1186785.1186792).
"""
function subdivide_simplex(s::Simplex{D,T,N}) where {D,T,N}
    color_schemes = subdivide_simplex_scheme(Val(D))
    return map(color_schemes) do color_scheme
        Simplex{D,T,N}(
            SVector{N}(
                (s.vertices[i] + s.vertices[j]) / 2 for (i, j) in eachcol(color_scheme)
            ),
        )
    end
end
