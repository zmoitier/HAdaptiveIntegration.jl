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

include("domain_subdivision.jl")

# Tabulated cubature for supported domains
include("quadrature_rules.jl")

# Subdivision strategies for various domains
include("quadrature_embedded.jl")

include("integrate.jl")
export integrate

include("quadrature_check.jl")

end
