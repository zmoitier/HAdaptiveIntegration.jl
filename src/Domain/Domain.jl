module Domain

# Domain types
export AbstractDomain, Simplex, Triangle, Tetrahedron, Orthotope, Rectangle, Cuboid
# Methods on domains
export dimension, reference_domain, map_from_reference, map_to_reference, abs_det_jac
# Subdivision methods
export subdivide_simplex,
    subdivide_triangle,
    subdivide_tetrahedron,
    subdivide_orthotope,
    subdivide_rectangle,
    subdivide_cuboid

using ..HAdaptiveIntegration: SVector, det, norm, promote_to_float, setindex

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
    dimension(::Type{<:AbstractDomain{D}}) where {D}

Return the dimension `D` of the given domain `DOM`.
"""
dimension(::Type{<:AbstractDomain{D}}) where {D} = D

# D-dimensional domains
include("simplex.jl")
include("orthotope.jl")

# 2-dimensional domains
include("triangle.jl")
include("rectangle.jl")

# 3-dimensional domains
include("tetrahedron.jl")
include("cuboid.jl")

"""
    reference_domain(::Type{<:AbstractDomain})

Return the reference domain for the given domain type. See [`reference_simplex`](@ref) and
[`reference_orthotope`](@ref) for more details.
"""
reference_domain(::Type{<:Simplex{D,T}}) where {D,T} = reference_simplex(D, T)
reference_domain(::Type{<:Simplex{D}}) where {D} = reference_simplex(D)

reference_domain(::Type{<:Orthotope{D,T}}) where {D,T} = reference_orthotope(D, T)
reference_domain(::Type{<:Orthotope{D}}) where {D} = reference_orthotope(D)

end
