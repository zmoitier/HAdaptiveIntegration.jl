module Rule

export AbstractRule, TabulatedEmbeddedCubature, orders, EmbeddedCubature, embedded_cubature
# Rules for a simplex
export GrundmannMoeller, RadonLaurie, TRIANGLE_GM19, TRIANGLE_RL19, TETRAHEDRON_GM35
# Rules for an orthotope
export GenzMalik, SQUARE_GM17, SQUARE_CH21, SQUARE_CH25, CUBE_GM33, CUBE_BE65, CUBE_BE115
# Rules for a segment
export SEGMENT_GK7, SEGMENT_GK15, SEGMENT_GK31

using ..HAdaptiveIntegration.Domain
using ..HAdaptiveIntegration: SVector, norm, promote_to_float

"""
    abstract type AbstractRule{DOM<:AbstractDomain}

Abstract type for a cubature rule on a domain `DOM`.

## Type Parameters:
- `DOM`: `[reference_domain](@ref)(DOM)` gives the reference domain on which the embedded
   cubature is assume to be set.

## Mandatory methods:
- [`embedded_cubature`](@ref)

## Useful (but non-mandatory) methods:
- [`orders`](@ref)
"""
abstract type AbstractRule{DOM<:AbstractDomain} end

include("tabulated.jl")
include("embedded_cubature.jl")

include("simplex.jl")
include("orthotope.jl")

include("segment.jl")
include("triangle.jl")
include("square.jl")
include("tetrahedron.jl")
include("cube.jl")

end
