module AdaptiveSimplexQuadrature

using DataStructures
using LinearAlgebra
using StaticArrays

# Supported integration domains
include("simplex.jl")
include("hyperrectangle.jl")

# Tabulated quadratures for supported domains
include("quadrature_rules.jl")

# Subdivision strategies for various domains
include("subdivision.jl")
include("quadrature_embedded.jl")

include("integrate.jl")

include("quadrature_check.jl")

end
