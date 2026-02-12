module Domain

# Domain types
export AbstractDomain, Segment, Triangle, Rectangle, Tetrahedron, Cuboid, Simplex, Orthotope
# Methods on domains
export dimension, reference_domain, map_from_reference, map_to_reference
# Subdivision methods
export subdivide_segment,
    subdivide_simplex,
    subdivide_triangle,
    subdivide_tetrahedron,
    subdivide_orthotope,
    subdivide_rectangle,
    subdivide_cuboid

using ..HAdaptiveIntegration.LinearAlgebra: det, norm
using ..HAdaptiveIntegration.StaticArrays: MMatrix, MVector, SMatrix, SVector, setindex

# AbstractDomain type
include("abstract_domain.jl")

# D-dimensional domains
include("simplex.jl")
include("orthotope.jl")

# 1-dimensional domains
include("segment.jl")

# 2-dimensional domains
include("triangle.jl")
include("rectangle.jl")

# 3-dimensional domains
include("tetrahedron.jl")
include("cuboid.jl")

reference_domain(::Type{<:Segment{T}}) where {T} = reference_segment(T)
reference_domain(::Type{<:Segment}) = reference_segment()

reference_domain(::Type{<:Simplex{D,T}}) where {D,T} = reference_simplex(Val(D), T)
reference_domain(::Type{<:Simplex{D}}) where {D} = reference_simplex(Val(D))

reference_domain(::Type{<:Orthotope{D,T}}) where {D,T} = reference_orthotope(Val(D), T)
reference_domain(::Type{<:Orthotope{D}}) where {D} = reference_orthotope(Val(D))

end
