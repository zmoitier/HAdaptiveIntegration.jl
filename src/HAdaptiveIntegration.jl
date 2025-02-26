module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

"""
    abstract type Domain{D,T}

Abstract type for integration domains in `D` dimensions with type `T<:Real`.
"""
abstract type Domain{D,T} end

"""
    map_from_reference(d::Domain)::Function

Return an anonymous function that maps the reference domain to the physical domain `d`.
"""
function map_from_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is unimplemented for type $(typeof(d))."
end

"""
    map_to_reference(d::Domain)::Function

Return an anonymous function that maps the physical domain `d` to the reference domain.
"""
function map_to_reference(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_to_reference` is unimplemented for type $(typeof(d))."
end

"""
    abs_det_jacobian(d::Domain)

The absolute value of the Jacobian's determinant of the map from the reference domain to the physical domain `d`.
"""
function abs_det_jacobian(d::Domain{D,T}) where {D,T<:Real}
    @error "`map_from_reference` is unimplemented for $(typeof(d))."
end

# Supported integration domains
include("domain_simplex.jl")
export triangle, reference_triangle, tetrahedron, reference_tetrahedron

include("domain_orthotope.jl")
export segment, reference_segment, rectangle, reference_rectangle, cuboid, reference_cuboid

# Subdivision strategies for various domains
include("domain_subdivision.jl")

const TABULATED_SUBDIVISION = Dict(
    :segment => [subdivide_segment2],
    :triangle => [subdivide_triangle2, subdivide_triangle4],
    :rectangle => [subdivide_rectangle4],
    :tetrahedron => [subdivide_tetrahedron8],
    :cuboid => [subdivide_cuboid8],
)

function default_subdivision(d::Domain)
    @error "no default subdivision for $(typeof(d))."
end

default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cuboid) = subdivide_cuboid8

# Embedded cubature
include("cubature_embedded.jl")
export EmbeddedCubatureRaw, embedded_cubature

include("cubature_check.jl")

# Tabulated cubature rule for supported domains
include("rule_segment.jl")
include("rule_triangle.jl")
include("rule_square.jl")
include("rule_tetrahedron.jl")
include("rule_cube.jl")

const TABULATED_EMBEDDED_CUBATURE = Dict(
    :segment => ["SEGMENT_G7K15", "SEGMENT_G15K31"],
    :triangle => ["TRIANGLE_LAURIE_RADON"],
    :rectangle => ["SQUARE_COOLS_HAEGEMANS"],
    :tetrahedron => [],
    :cuboid => [],
)

function default_embedded_cubature(d::Domain)
    @error "no default embedded cubature for $(typeof(d))."
end

function default_embedded_cubature(::Segment{T}) where {T}
    return embedded_cubature_from_raw(SEGMENT_G7K15, T)
end
function default_embedded_cubature(::Triangle{T}) where {T}
    return embedded_cubature_from_raw(TRIANGLE_LAURIE_RADON, T)
end
function default_embedded_cubature(::Rectangle{T}) where {T}
    return embedded_cubature_from_raw(SQUARE_COOLS_HAEGEMANS, T)
end
# default_embedded_cubature(::Tetrahedron) = 
# default_embedded_cubature(::Cuboid) = 

# Compute integrals
include("integrate.jl")
export integrate

end
