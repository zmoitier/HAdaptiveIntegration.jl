module AdaptiveSimplexQuadrature

import Base: length, size

using DataStructures
using LinearAlgebra
using StaticArrays

include("simplex.jl")

include("quadrature.jl")
include("quadrature_embedded.jl")
include("quadrature_rules.jl")
include("quadrature_check.jl")

include("integrate.jl")

end
