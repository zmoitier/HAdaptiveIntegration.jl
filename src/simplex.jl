"""
    struct Simplex{N,T,Np1}

A simplex in `N` dimensions with Np1=N+1 points of type `T`.
"""
struct Simplex{N,T,Np1}
    points::SVector{Np1,SVector{N,T}}
end

function Simplex(pts...)
    N = length(pts) - 1 # ambient dimension
    @assert all(pt -> length(pt) == N, pts)
    pts_vec = SVector(SVector.(pts))
    return Simplex(pts_vec)
end

function Simplex{N,T,Np1}(pts...)::Simplex{N,T,Np1} where {N,T,Np1}
    return Simplex(pts...)
end

const Segment{T}     = Simplex{1,T,2}
const Triangle{T}    = Simplex{2,T,3}
const Tetrahedron{T} = Simplex{3,T,4}

# default types
Segment(args...) = Segment{Float64}(args...)
Triangle(args...) = Triangle{Float64}(args...)
Tetrahedron(args...) = Tetrahedron{Float64}(args...)

"""
    reference_simplex(N::Int, T::DataType)

Return the reference N-simplex. The reference simplex has vertices given by
`(0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1)`.
"""
function reference_simplex(N::Int, T::DataType)
    @assert N ≥ 1 "N = $N must be greater than 1."

    vertices = [zeros(T, N)]
    for i in 1:N
        tmp = zeros(T, N)
        tmp[i] = T(1)
        push!(vertices, tmp)
    end

    Np1 = N + 1
    return Simplex{N,T,Np1}(SVector{Np1}(SVector{N}.(vertices)))
end

"""
    map_from_ref(s::Simplex)

Return an anonymous function that maps the reference simplex to the physical simplex.
The reference simplex has vertices given by
`(0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1)`.
"""
function map_from_ref(s::Simplex{N,T,Np1})::Function where {N,T,Np1}
    return u ->
        (1 - sum(u)) * s.points[1] + sum(c * p for (c, p) in zip(u, s.points[2:Np1]))
end

"""
    map_to_ref(s::Simplex)

Return an anonymous function that maps the physical simplex to the reference simplex.
The reference simplex has vertices given by
`(0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1)`.
"""
function map_to_ref(s::Simplex{N,T,Np1})::Function where {N,T,Np1}
    e0 = s.points[1]

    M = zeros(T, (N, N))
    for i in 1:N
        M[:, i] .= s.points[i+1] - e0
    end

    M_inv = SMatrix{N,N}(inv(M))
    return u -> M_inv * (u - e0)
end

"""
    measure(s::Simplex)

The n-dimensional volume of the simplex.
"""
function measure(s::Simplex{N,T,Np1}) where {N,T,Np1}
    v = s.points
    return abs(det(reinterpret(reshape, T, [x - v[1] for x in v[2:Np1]]))) / factorial(N)
end

"""
    measure(s::Segment)

The length of the segment.
"""
function measure(s::Segment{T}) where {T}
    return abs(s.points[1][1] - s.points[2][1])
end

"""
    measure(t::Triangle)

The area of the triangle.
"""
function measure(t::Triangle{T}) where {T}
    v = t.points
    e1 = (v[2] - v[1])
    e2 = (v[3] - v[1])
    return norm(cross(e1, e2)) / 2
end

"""
    measure(t::Tetrahedron)

The volume of the Tetrahedron.
"""
function measure(t::Tetrahedron{T}) where {T}
    v = t.points
    e1 = (v[2] - v[1])
    e2 = (v[3] - v[1])
    e3 = (v[4] - v[1])
    return norm(dot(e1, cross(e2, e3))) / 6
end

function subdivide(s::Segment{T}) where {T}
    a, b = s.points
    m = (a + b) / 2
    return (Segment(a, m), Segment(m, b))
end

function subdivide(t::Triangle{T}) where {T}
    p1, p2, p3 = t.points
    p12 = (p1 + p2) / 2
    p23 = (p2 + p3) / 2
    p31 = (p3 + p1) / 2
    return (
        Triangle(p1, p12, p31),
        Triangle(p2, p23, p12),
        Triangle(p3, p31, p23),
        Triangle(p12, p23, p31),
    )
end
