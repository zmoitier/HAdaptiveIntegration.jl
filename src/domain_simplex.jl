"""
    struct Simplex{D,T,N}

A simplex in `D` dimensions with N=D+1 points of type `T`.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    points::SVector{N,SVector{D,T}}
end

"""
    simplex(points...)

A simplex in `D` dimensions with N=D+1 points of type `T`.
"""
function simplex(points...)
    D = length(points) - 1
    @assert all(pt -> length(pt) == D, points)
    N = D + 1

    vec_pts = SVector{N}(SVector{D}.(points))
    return Simplex(vec_pts)
end

"""
    reference_simplex(D, T=Float64)

Return the reference D-simplex given by the vertices `(0,...,0), (1,0,...,0), (0,1,0,...,0), (0,...,0,1)`.
"""
function reference_simplex(D::Int, T::DataType=Float64)
    @assert D ≥ 1 "D = $D must be greater than 1."
    N = D + 1

    points = [zeros(T, D) for _ in 1:N]
    for i in 1:D
        points[i + 1][i] = 1
    end

    return Simplex(SVector{N,SVector{D,T}}(SVector{D,T}.(points)))
end

function map_from_reference(s::Simplex{D,T,N}) where {D,T<:Real,N}
    return u -> (1 - sum(u)) * s.points[1] + sum(c * p for (c, p) in zip(u, s.points[2:N]))
end

function map_to_reference(s::Simplex{D,T,N}) where {D,T<:Real,N}
    v = s.points
    M = inv(SMatrix{D,D}(reinterpret(reshape, T, [x - v[1] for x in v[2:N]])))
    return u -> M * (u - v[1])
end

function abs_det_jacobian(s::Simplex{D,T,N}) where {D,T<:Real,N}
    v = s.points
    return abs(det(SMatrix{D,D}(reinterpret(reshape, T, [x - v[1] for x in v[2:N]]))))
end

"""
A triangle in 2 dimensions with 3 vertices of type `T`.
"""
const Triangle{T} = Simplex{2,T,3}

function Triangle(a::SVector{2,T}, b::SVector{2,T}, c::SVector{2,T}) where {T<:Real}
    return Simplex(SVector{3,SVector{2,T}}([a, b, c]))
end

"""
    triangle(a, b, c)

A triangle in 2 dimensions given by the 2d-points `a`, `b`, and `c`.
"""
function triangle(a, b, c)
    return simplex(a, b, c)
end

"""
    reference_triangle(T=Float64)

Return the reference triangle given by the points `(0,0), (1,0), (0,1)`.
"""
function reference_triangle(T::DataType=Float64)
    return reference_simplex(2, T)
end

function abs_det_jacobian(t::Triangle{T}) where {T<:Real}
    e1 = t.points[2] - t.points[1]
    e2 = t.points[3] - t.points[1]
    return norm(cross(e1, e2))
end

"""
A tetrahedron in 3 dimensions with 4 points of type `T`.
"""
const Tetrahedron{T} = Simplex{3,T,4}

function Tetrahedron(
    a::SVector{3,T}, b::SVector{3,T}, c::SVector{3,T}, d::SVector{3,T}
) where {T<:Real}
    return Simplex(SVector{4,SVector{3,T}}([a, b, c, d]))
end

"""
    tetrahedron(a, b, c, d)

A tetrahedron in 3 dimensions given by the 3d-points `a`, `b`, `c`, and `d`.
"""
function tetrahedron(a, b, c, d)
    return simplex(a, b, c, d)
end

"""
    reference_tetrahedron(T::DataType=Float64)

Return the reference tetrahedron given by the vertices `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`.
"""
function reference_tetrahedron(T::DataType=Float64)
    return reference_simplex(3, T)
end

function abs_det_jacobian(t::Tetrahedron{T}) where {T<:Real}
    e1 = t.points[2] - t.points[1]
    e2 = t.points[3] - t.points[1]
    e3 = t.points[4] - t.points[1]
    return abs(dot(e1, cross(e2, e3)))
end
