module Domain

# Domain types
export AbstractDomain,
    Simplex,
    Triangle,
    Tetrahedron,
    Orthotope,
    Segment,
    Rectangle,
    Cuboid,
    simplex,
    triangle,
    tetrahedron,
    orthotope,
    segment,
    rectangle,
    cuboid
# Methods on domains
export dimension, reference_domain, map_from_reference, map_to_reference, abs_det_jac
# Subdivision methods
export subdivide_simplex,
    subdivide_triangle2,
    subdivide_triangle4,
    subdivide_tetrahedron8,
    subdivide_orthotope,
    subdivide_segment2,
    subdivide_rectangle4,
    subdivide_cuboid8

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
    dimension(::Type{DOM}) where {DOM<:AbstractDomain{D}}

Return the dimension `D` of the given domain `DOM`.
"""
dimension(::Type{AbstractDomain{D,T}}) where {D,T} = D
dimension(::Type{AbstractDomain{D}}) where {D} = D
dimension(::Type{DOM}) where {DOM<:AbstractDomain} = dimension(supertype(DOM))

include("simplex.jl")
include("orthotope.jl")

"""
    reference_domain(::Type{<:AbstractDomain})

Return the reference domain for the given domain type.
"""
reference_domain(::Type{Simplex{D,N,T}}) where {D,N,T} = reference_simplex(T, D)
reference_domain(::Type{Simplex{D,N}}) where {D,N} = reference_simplex(float(Int), D)
reference_domain(::Type{Simplex{D}}) where {D} = reference_simplex(float(Int), D)

reference_domain(::Type{Orthotope{D,T}}) where {D,T} = reference_orthotope(T, D)
reference_domain(::Type{Orthotope{D}}) where {D} = reference_orthotope(float(Int), D)

end
