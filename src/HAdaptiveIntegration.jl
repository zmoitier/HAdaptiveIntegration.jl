module HAdaptiveIntegration

using DataStructures
using LinearAlgebra
using StaticArrays

using GrundmannMoeller: grundmann_moeller

# Supported integration domains
include("domain.jl")
export Orthotope, segment, rectangle, cuboid, simplex, triangle, tetrahedron

# Subdivision strategies for various domains
include("domain_subdivision.jl")

# Embedded cubature
include("cubature_embedded.jl")
include("cubature_check.jl")

# Tabulated cubature rule for supported domains
include("rule_segment.jl")
include("rule_triangle.jl")
include("rule_square.jl")
include("rule_tetrahedron.jl")
include("rule_cube.jl")

const LIST_EMBEDDED_CUBATURE = [
    "segment" => ["SEGMENT_G7K15", "SEGMENT_G15K31"],
    "rectangle" => ["SQUARE_CH21G25"],
    "triangle" => ["TRIANGLE_R7L19"],
    "cuboid" => ["CUBE_BE65"],
    "tetrahedron" => [],
]

function default_embedded_cubature(d::Domain)
    @error "no default embedded cubature for $(typeof(d))."
end
@generated function default_embedded_cubature(::Segment{T}) where {T}
    ec = embedded_cubature(SEGMENT_G7K15, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Rectangle{T}) where {T}
    ec = embedded_cubature(SQUARE_CH21G25, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Triangle{T}) where {T}
    ec = embedded_cubature(TRIANGLE_R7L19, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Cuboid{T}) where {T}
    ec = embedded_cubature(CUBE_BE65, T)
    return :($ec)
end
@generated function default_embedded_cubature(::Simplex{D,T,N}) where {D,T,N}
    ec = embedded_cubature(GrundmannMoeller, D, 7, T)
    return :($ec)
end

# Compute integrals
include("integrate.jl")
export integrate

end
