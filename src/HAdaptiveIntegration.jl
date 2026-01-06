module HAdaptiveIntegration

using DataStructures: BinaryHeap
using LinearAlgebra: LinearAlgebra
using StaticArrays: StaticArrays

# Supported integration domains
include("Domain/Domain.jl")
using .Domain
export Segment, Triangle, Rectangle, Tetrahedron, Cuboid, Simplex, Orthotope

# Tabulated cubature rule for supported domains
include("Rule/Rule.jl")
using .Rule

# Default subdivision and embedded cubature for supported domain
include("default.jl")

# Compute integrals
include("integrate.jl")
export integrate

# Extensions
include("extension.jl")

end
