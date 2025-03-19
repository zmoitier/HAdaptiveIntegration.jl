"""
    abstract type AbstractRule{DOM<:AbstractDomain}

Abstract type for a cubature rule on a domain `DOM`.

## Mandatory methods:
- [`embedded_cubature`](@ref)
"""
abstract type AbstractRule{DOM<:AbstractDomain} end

"""
    @kwdef struct TabulatedEmbeddedCubature

An embedded cubature rule consisting of a high order cubature rule and a low order cubature
rule. Note that the low order cubature uses `nodes[1:L]` as its nodes where `L` is the
length of the `weights_low`. The cubature nodes and weights are assume to be for the
reference domain (use the [`reference_domain`](@ref) function to get the reference domain).

## Fields:
- `description::String`: description of the embedded cubature.
- `reference::String`: where the values are found.
- `nb_significant_digits::Int`: number of significant digits on the node and weight values,
  `10^-nb_significant_digits` give the relative precision of the values.
- `nodes::Vector{Vector{String}}`: the cubature nodes.
- `weights_high::Vector{String}`: the cubature weights for the high order cubature.
- `order_high::Int`: order of the high order cubature.
- `weights_low::Vector{String}`: the cubature weights for the low order cubature.
- `order_low::Int`: order of the low order cubature.

## Invariants (check at construction):
- `length(nodes) == length(weights_high)`
- `length(weights_high) ≥ length(weights_low)`
- `order_high ≥ order_low`
"""
struct TabulatedEmbeddedCubature{D} # {D,DOM<:AbstractDomain{D,T}} <: AbstractRule{DOM}
    description::String
    reference::String
    nb_significant_digits::Int
    nodes::Vector{NTuple{D,String}}
    weights_high::Vector{String}
    order_high::Int
    weights_low::Vector{String}
    order_low::Int

    function TabulatedEmbeddedCubature{D}(;
        description::String,
        reference::String,
        nb_significant_digits::Int,
        nodes::Vector{NTuple{D,String}},
        weights_high::Vector{String},
        order_high::Int,
        weights_low::Vector{String},
        order_low::Int,
    ) where {D} # {DOM<:AbstractDomain}
        @assert length(nodes) == length(weights_high)
        @assert length(weights_high) ≥ length(weights_low)
        @assert order_high ≥ order_low
        return new{D}(
            description,
            reference,
            nb_significant_digits,
            nodes,
            weights_high,
            order_high,
            weights_low,
            order_low,
        )
    end
end

"""
   struct GrundmannMoeller

Cubature rule for a `dim`-simplex of degree `deg`.
"""
struct GrundmannMoeller # {D,N,T} <: AbstractRule{Simplex{D,N,T}}
    dim::Int
    deg::Int
end
