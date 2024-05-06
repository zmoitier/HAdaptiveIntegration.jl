module AdaptiveSimplexQuadrature

using StaticArrays

include("quadrules.jl")
include("simplex.jl")
include("embeddedquadrature.jl")
include("adaptive.jl")

end
