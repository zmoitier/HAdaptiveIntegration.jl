"""
    struct Simplex{N,T,Np1}

A simplex in `N` dimensions with points of type `T`.
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

function (t::Triangle)(u)
    v = t.points
    return v[1] + (v[2] - v[1]) * u[1] + (v[3] - v[1]) * u[2]
end

function subdivide(t::Triangle{T}) where {T}
    p1, p2, p3 = t.points
    p12 = (p1 + p2) / 2
    p23 = (p2 + p3) / 2
    p31 = (p3 + p1) / 2
    return (
        Triangle{T}(p1, p12, p31),
        Triangle{T}(p2, p23, p12),
        Triangle{T}(p3, p31, p23),
        Triangle{T}(p12, p23, p31),
    )
end
