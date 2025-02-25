module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

"""
    abstract type Domain{D,T}

Abstract type for integration domains in `D` dimensions.
"""
abstract type Domain{D,T} end

# Supported integration domains
include("domain_simplex.jl")
export triangle, reference_triangle, tetrahedron, reference_tetrahedron

include("domain_orthotope.jl")
export segment, reference_segment, rectangle, reference_rectangle, cuboid, reference_cuboid

# Subdivision strategies for various domains
include("domain_subdivision.jl")

function default_subdivision(domain::Domain)
    @error "no default subdivision for $(typeof(domain))."
end

default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cuboid) = subdivide_cuboid8

# Embedded cubature
include("cubature_embedded.jl")
export embedded_cubature

include("cubature_check.jl")

# Tabulated cubature rule for supported domains
include("rule_segment.jl")
include("rule_triangle.jl")
include("rule_square.jl")
include("rule_tetrahedron.jl")
include("rule_cube.jl")

const TABULATED_EMBEDDED_CUBATURE = Dict(
    :segment => [SEGMENT_G7K15, SEGMENT_G15K31],
    :triangle => [],
    :square => [],
    :tetrahedron => [],
    :cube => [],
)

function default_embedded_cubature(domain::Domain)
    @error "no default embedded cubature for $(typeof(domain))."
end

function default_embedded_cubature(::Segment{T}) where {T}
    return embedded_cubature_from_raw(SEGMENT_G7K15, T)
end
# default_embedded_cubature(::Triangle) = 
# default_embedded_cubature(::Rectangle) = 
# default_embedded_cubature(::Tetrahedron) = 
# default_embedded_cubature(::Cuboid) = 

# include("integrate.jl")
# export integrate

end
