module Rule

export AbstractRule, TabulatedEmbeddedCubature, orders, EmbeddedCubature, embedded_cubature
# Rules for a segment
export SEGMENT_GK7, SEGMENT_GK15, SEGMENT_GK31
# Rules for a simplex
export GrundmannMoeller, RadonLaurie, TRIANGLE_GM19, TRIANGLE_RL19, TETRAHEDRON_GM35
# Rules for an orthotope
export GenzMalik, SQUARE_GM17, SQUARE_CH21, SQUARE_CH25, CUBE_GM33, CUBE_BE65, CUBE_BE115

using ..HAdaptiveIntegration.Domain
using LinearAlgebra: norm
using StaticArrays: SVector

"""
    abstract type AbstractRule{DOM<:AbstractDomain}

Abstract type for a cubature rule on a domain `DOM`.

## Type Parameters:
- `DOM`: [`reference_domain(DOM)`](@ref) gives the reference domain on which the embedded
   cubature is assume to be set.

## Mandatory methods:
- [`embedded_cubature`](@ref)

## Useful (but non-mandatory) methods:
- [`orders`](@ref)
"""
abstract type AbstractRule{DOM<:AbstractDomain} end

include("embedded_cubature.jl")
include("tabulated.jl")

# D-dimensional rules
include("simplex.jl")
include("orthotope.jl")

# 1-dimensional rules
include("segment.jl")

# 2-dimensional rules
include("triangle.jl")
include("square.jl")

# 3-dimensional rules
include("tetrahedron.jl")
include("cube.jl")

end
