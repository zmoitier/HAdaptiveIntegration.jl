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
export Triangle, Tetrahedron, reference_triangle, reference_tetrahedron

include("domain_orthotope.jl")
export Segment, Rectangle, Cuboid, reference_segment, reference_rectangle, reference_cuboid

include("domain_subdivision.jl")

# Tabulated cubature for supported domains
include("quadrature_rules.jl")

# Subdivision strategies for various domains
include("quadrature_embedded.jl")

include("integrate.jl")
export integrate

include("quadrature_check.jl")

end
