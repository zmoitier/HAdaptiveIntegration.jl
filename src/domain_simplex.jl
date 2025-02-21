"""
    struct Simplex{D,T,N} <: Domain{D,T}

A simplex in `D` dimensions with N=D+1 points of type `T`.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    points::SVector{N,SVector{D,T}}
end

function Simplex(points...)
    D = length(points) - 1
    @assert all(pt -> length(pt) == D, points)
    N = D + 1

    vec_pts = SVector{N}(SVector{D}.(points))
    return Simplex(vec_pts)
end

function Simplex{D,T,N}(points...) where {D,T,N}
    @assert all(pt -> length(pt) == D, points)

    vec_pts = SVector(SVector{D,T}.(points))
    return Simplex(vec_pts)
end

"""
    reference_simplex(D::Int, T::DataType)

Return the reference D-simplex given by the points `(0,...,0), (1,0,...,0), (0,1,0,...,0),
(0,...,0,1)`.
"""
function reference_simplex(D::Int, T::DataType=Float64)
    @assert D ≥ 1 "D = $D must be greater than 1."
    N = D + 1

    points = [zeros(T, D) for _ in 1:N]
    for i in 1:D
        points[i + 1][i] = 1
    end

    return Simplex{D,T,N}(SVector{N}(SVector{D}.(points)))
end

"""
    map_from_reference(s::Simplex)::Function

Return an anonymous function that maps the reference simplex to the physical simplex `s`.
"""
function map_from_reference(s::Simplex{D,T,N})::Function where {D,T,N}
    return u -> (1 - sum(u)) * s.points[1] + sum(c * p for (c, p) in zip(u, s.points[2:N]))
end

"""
    abs_jacobian_determinant(s::Simplex)

The determinant of the Jacobian of the map from the reference simplex to the physical simplex `s`.
"""
function abs_jacobian_determinant(s::Simplex{D,T,N}) where {D,T,N}
    v = s.points
    return abs(det(SMatrix{D,D}(reinterpret(reshape, T, [x - v[1] for x in v[2:N]]))))
end

# Special simplex
const Triangle{T} = Simplex{2,T,3}
const Tetrahedron{T} = Simplex{3,T,4}

# default types
Triangle(points...) = Triangle{Float64}(points...)
Tetrahedron(points...) = Tetrahedron{Float64}(points...)

"""
    reference_triangle(T::DataType=Float64)

Return the reference triangle given by the points `(0,0), (1,0), (0,1)`.
"""
function reference_triangle(T::DataType=Float64)
    return reference_simplex(2, T)
end

"""
    reference_tetrahedron(T::DataType=Float64)

Return the reference tetrahedron given by the points `(0,0,0), (1,0,0), (0,1,0), (0,0,1)`.
"""
function reference_tetrahedron(T::DataType=Float64)
    return reference_simplex(3, T)
end

"""
    abs_jacobian_determinant(t::Triangle)

The determinant of the Jacobian of the map from the reference triangle to the physical triangle `t`.
"""
function abs_jacobian_determinant(t::Triangle{T}) where {T}
    e1 = t.points[2] - t.points[1]
    e2 = t.points[3] - t.points[1]
    return norm(cross(e1, e2))
end

"""
    abs_jacobian_determinant(t::Tetrahedron)

The determinant of the Jacobian of the map from the reference tetrahedron to the physical
tetrahedron `t`.
"""
function abs_jacobian_determinant(t::Tetrahedron{T}) where {T}
    e1 = t.points[2] - t.points[1]
    e2 = t.points[3] - t.points[1]
    e3 = t.points[4] - t.points[1]
    return norm(dot(e1, cross(e2, e3)))
end
