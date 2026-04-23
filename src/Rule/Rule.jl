module Rule

export AbstractRule, TabulatedEmbeddedCubature, orders, EmbeddedCubature, embedded_cubature
# Rules for a segment
export SEGMENT_GK7, SEGMENT_GK15, SEGMENT_GK31
# Rules for a simplex
export GrundmannMoeller, RadonLaurie
# Rules for an orthotope
export GenzMalik, SQUARE_CH21, SQUARE_CH25, CUBE_BE65, CUBE_BE115

using ..HAdaptiveIntegration.Domain
using ..HAdaptiveIntegration.LinearAlgebra: LinearAlgebra
using ..HAdaptiveIntegration.StaticArrays: SVector

"""
    abstract type AbstractRule{DOM<:AbstractDomain}

Abstract type for a cubature rule on a domain `DOM`.

## Type Parameters:
- `DOM`: [`reference_domain(DOM)`](@ref) gives the reference domain on which the embedded
  cubature is assumed to be defined.

## Mandatory methods:
- [`embedded_cubature`](@ref)

## Useful (but non-mandatory) methods:
- [`orders`](@ref)
"""
abstract type AbstractRule{DOM <: AbstractDomain} end

"""
    embedded_cubature(ar::AbstractRule, T=float(Int))

    embedded_cubature(nodes, weights_high, weights_low, T=float(Int))

Construct an embedded cubature with element type `T`.

The constructor can be called from a subtype of [`AbstractRule`](@ref), or from explicit
`nodes`, `weights_high`, and `weights_low` data. Available rule types include:
- [`TabulatedEmbeddedCubature`](@ref)
- [`RadonLaurie`](@ref)
- [`GrundmannMoeller`](@ref)
- [`GenzMalik`](@ref)
"""
function embedded_cubature end

"""
    orders(rule::AR) where {AR<:AbstractRule}

Return `(order_high, order_low)` for `rule`.
"""
function orders end

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
