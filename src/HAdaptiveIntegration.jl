module HAdaptiveIntegration

using DataStructures: BinaryHeap
using LinearAlgebra: det, norm
using StaticArrays: SVector, setindex

# Supported integration domains
include("Domain/Domain.jl")
using .Domain
export Simplex, Triangle, Tetrahedron, Orthotope, Rectangle, Cuboid

# Tabulated cubature rule for supported domains
include("Rule/Rule.jl")
using .Rule

# Default subdivision and embedded cubature for supported domain
include("default.jl")

# Compute integrals
include("integrate.jl")
export integrate

end
