module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

# Abstract domain type
include("domain.jl")
export reference_domain

# Supported integration domains
include("domain_simplex.jl")
export Triangle, triangle, Tetrahedron, tetrahedron, Simplex, simplex, reference_simplex

include("domain_orthotope.jl")
export Segment,
    segment, Rectangle, rectangle, Cuboid, cuboid, Orthotope, orthotope, reference_orthotope

const LIST_DOMAIN_TYPE = [
    "dimension 1" => ["segment"],
    "dimension 2" => ["rectangle", "triangle"],
    "dimension 3" => ["cuboid", "tetrahedron"],
    "dimension D" => ["orthotope", "simplex"],
]

# Subdivision strategies for various domains
include("domain_subdivision.jl")
export subdivide_segment2,
    subdivide_segment3,
    subdivide_triangle2,
    subdivide_triangle4,
    subdivide_rectangle4,
    subdivide_tetrahedron8,
    subdivide_cuboid8,
    check_subdivision

const LIST_SUBDIVISION_ALGO = [
    "segment" => ["subdivide_segment2", "subdivide_segment3"],
    "rectangle" => ["subdivide_rectangle4"],
    "triangle" => ["subdivide_triangle2", "subdivide_triangle4"],
    "cuboid" => ["subdivide_cuboid8"],
    "tetrahedron" => ["subdivide_tetrahedron8"],
]

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
export EmbeddedCubatureRaw, embedded_cubature_from_raw, embedded_cubature

include("cubature_check.jl")
export check_order

# Tabulated cubature rule for supported domains
include("rule_segment.jl")
export SEGMENT_G7K15, SEGMENT_G15K31

include("rule_triangle.jl")
export TRIANGLE_R7_L19

include("rule_square.jl")
export SQUARE_CH21_G25

include("rule_tetrahedron.jl")
# export

include("rule_cube.jl")
#export

const LIST_EMBEDDED_CUBATURE = [
    "segment" => ["SEGMENT_G7K15", "SEGMENT_G15K31"],
    "rectangle" => ["SQUARE_CH21_G25"],
    "triangle" => ["TRIANGLE_R7_L19"],
    "cuboid" => [],
    "tetrahedron" => [],
]

function default_embedded_cubature(d::Domain)
    @error "no default embedded cubature for $(typeof(d))."
end

function default_embedded_cubature(::Segment{T}) where {T}
    return embedded_cubature_from_raw(SEGMENT_G7K15, T)
end
function default_embedded_cubature(::Triangle{T}) where {T}
    return embedded_cubature_from_raw(TRIANGLE_R7_L19, T)
end
function default_embedded_cubature(::Rectangle{T}) where {T}
    return embedded_cubature_from_raw(SQUARE_CH21_G25, T)
end
# default_embedded_cubature(::Tetrahedron) = 
# default_embedded_cubature(::Cuboid) = 

# Compute integrals
include("integrate.jl")
export integrate

end
