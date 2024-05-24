module AdaptiveSimplexQuadrature

using DataStructures
using LinearAlgebra
using StaticArrays

include("quadrules.jl")
include("simplex.jl")
include("embeddedquadrature.jl")
include("integrate.jl")
include("reference_values.jl")

end
