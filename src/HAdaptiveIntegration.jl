module HAdaptiveIntegration

using DataStructures
using GrundmannMoeller: grundmann_moeller
using LinearAlgebra
using StaticArrays

# Supported integration domains
include("domain.jl")
export orthotope, segment, rectangle, cuboid, simplex, triangle, tetrahedron

# Subdivision strategies for various domains
include("domain_subdivision.jl")

# Tabulated cubature rule for supported domains
include("rule.jl")
include("rule_orthotope.jl")
include("rule_simplex.jl")

# Embedded cubature
include("cubature_embedded.jl")
include("cubature_check.jl")

# Default subdivision and embedded cubature for supported domain
include("default.jl")

# Compute integrals
include("integrate.jl")
export integrate

end
