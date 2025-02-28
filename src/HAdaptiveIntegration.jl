module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

# Supported integration domains
include("domain.jl")
export Orthotope,
    Segment,
    segment,
    Rectangle,
    rectangle,
    Cuboid,
    cuboid,
    Triangle,
    triangle,
    Tetrahedron,
    tetrahedron,
    Simplex,
    simplex,
    LIST_DOMAIN_TYPE,
    reference_domain,
    map_from_reference,
    map_to_reference

# Subdivision strategies for various domains
include("domain_subdivision.jl")
export check_subdivision,
    subdivide_segment2,
    subdivide_segment3,
    subdivide_triangle2,
    subdivide_triangle4,
    subdivide_rectangle4,
    subdivide_tetrahedron8,
    subdivide_cuboid8,
    LIST_SUBDIVISION_ALGO,
    default_subdivision

# Embedded cubature
include("cubature_embedded.jl")
export EmbeddedCubatureRaw, embedded_cubature_from_raw, embedded_cubature

include("cubature_check.jl")
export check_order

# Tabulated cubature rule for supported domains
include("rule_segment.jl")
export SEGMENT_G7K15, SEGMENT_G15K31

include("rule_triangle.jl")
export TRIANGLE_R7L19

include("rule_square.jl")
export SQUARE_CH21G25

include("rule_tetrahedron.jl")
# export

include("rule_cube.jl")
#export

const LIST_EMBEDDED_CUBATURE = [
    "segment" => ["SEGMENT_G7K15", "SEGMENT_G15K31"],
    "rectangle" => ["SQUARE_CH21G25"],
    "triangle" => ["TRIANGLE_R7L19"],
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
    return embedded_cubature_from_raw(TRIANGLE_R7L19, T)
end
function default_embedded_cubature(::Rectangle{T}) where {T}
    return embedded_cubature_from_raw(SQUARE_CH21G25, T)
end
# default_embedded_cubature(::Tetrahedron) = 
# default_embedded_cubature(::Cuboid) = 

# Compute integrals
include("integrate.jl")
export integrate

end
