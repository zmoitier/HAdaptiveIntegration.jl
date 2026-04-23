module Domain

# Domain types
export AbstractDomain, Segment, Triangle, Rectangle, Tetrahedron, Cuboid, Simplex, Orthotope
# Methods on domains
export reference_domain, map_from_reference, map_to_reference
# Subdivision methods
export subdivide_segment,
    subdivide_simplex,
    subdivide_triangle,
    subdivide_tetrahedron,
    subdivide_orthotope,
    subdivide_rectangle,
    subdivide_cuboid

using ..HAdaptiveIntegration.LinearAlgebra: det
using ..HAdaptiveIntegration.StaticArrays: MMatrix, MVector, SMatrix, SVector, setindex

"""
    abstract type AbstractDomain{D,T}

Abstract type for integration domains in `D` dimensions with element type `T`.

## Type Parameters:
- `D`: dimension of the domain.
- `T`: element type of the domain.

## Mandatory methods:
- [`map_from_reference`](@ref)

## Useful (but non-mandatory) methods:
- [`reference_domain`](@ref)
- [`map_to_reference`](@ref)
"""
abstract type AbstractDomain{D, T} end

"""
    Φ, μ = map_from_reference(domain::DOM) where {DOM<:AbstractDomain}

Return `(Φ, μ)`, where `Φ` is an anonymous function that maps the reference domain to the
physical domain `domain`, and `μ` is the absolute value of the Jacobian determinant of `Φ`.
"""
function map_from_reference end

"""
    Ψ = map_to_reference(domain::DOM) where {DOM<:AbstractDomain}

Return `Ψ` an anonymous function that maps the physical `domain` to the reference domain.

## Constraints:
- For `Segment`, must have `xmax ≠ xmin`.
- For `Simplex{D}`, the vertices must form a valid `D`-dimensional simplex with non-zero
  volume.
- For `Orthotope`, must have `high_corner ≠ low_corner`.
"""
function map_to_reference end

"""
    reference_domain(::Type{<:AbstractDomain})

Return the reference domain for the given domain type. See [`reference_segment`](@ref),
[`reference_simplex`](@ref), and [`reference_orthotope`](@ref) for more details.
"""
function reference_domain end

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

reference_domain(::Type{<:Simplex{D, T}}) where {D, T} = reference_simplex(Val(D), T)
reference_domain(::Type{<:Simplex{D}}) where {D} = reference_simplex(Val(D))

reference_domain(::Type{<:Orthotope{D, T}}) where {D, T} = reference_orthotope(Val(D), T)
reference_domain(::Type{<:Orthotope{D}}) where {D} = reference_orthotope(Val(D))

end
