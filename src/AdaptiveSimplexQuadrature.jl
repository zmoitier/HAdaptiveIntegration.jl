module AdaptiveSimplexQuadrature

using DataStructures
using LinearAlgebra
using StaticArrays

"""
    abstract type Domain{N,T}

Abstract type for integration domains in `N` dimensions.
"""
abstract type Domain{N,T} end

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
