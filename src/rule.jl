"""
    abstract type AbstractRule{DOM<:AbstractDomain}

Abstract type for a cubature rule on a domain `DOM`.

## Type Parameters:
- `DOM`: `[reference_domain](@ref)(DOM)` gives the reference domain on which the embedded
   cubature is assume to be set.

## Mandatory methods:
- [`embedded_cubature`](@ref)
"""
abstract type AbstractRule{DOM<:AbstractDomain} end

"""
    @kwdef struct TabulatedEmbeddedCubature{DOM<:AbstractDomain} <: AbstractRule{DOM}

An embedded cubature rule consisting of a high order cubature rule and a low order cubature
rule. Note that the low order cubature uses `nodes[1:L]` as its nodes where `L` is the
length of the `weights_low`. The cubature nodes and weights are assumed to be for the
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
- `order_high ≥ order_low` ≥ 0
- `nb_significant_digits ≥ 0`
"""
struct TabulatedEmbeddedCubature{DOM<:AbstractDomain} <: AbstractRule{DOM}
    description::String
    reference::String
    nb_significant_digits::Int
    nodes::Vector{Vector{String}}
    weights_high::Vector{String}
    order_high::Int
    weights_low::Vector{String}
    order_low::Int

    function TabulatedEmbeddedCubature{DOM}(;
        description::String,
        reference::String,
        nb_significant_digits::Int,
        nodes::Vector{Vector{String}},
        weights_high::Vector{String},
        order_high::Int,
        weights_low::Vector{String},
        order_low::Int,
    ) where {DOM<:AbstractDomain}
        D = dimension(DOM)
        @assert all(n -> length(n) == D, nodes) "Each node must have length equal to the dimension D"
        @assert length(nodes) == length(weights_high) "The number of nodes must match the number of high-order weights"
        @assert length(weights_high) ≥ length(weights_low) "The length of high order weights must be greater or equal to the length of low-order weights"
        @assert order_high ≥ order_low ≥ 0 "order_high must be greater than or equal to order_low and orders must be non-negative"
        @assert nb_significant_digits ≥ 0
        return new{DOM}(
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
   struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}

Embedded cubature rule for a `D`-simplex.

## Type Parameters:
- `D`: The dimension of the simplex.

## Fields:
- `order_high::Int`: the high order of the rule.
- `order_low::Int`: the low order of the rule.

## Invariants (check at construction):
- `order_high` and `order_low` must be odd.
- must have `order_high > order_low ≥ 1`.
"""
struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}
    order_high::Int
    order_low::Int

    function GrundmannMoeller{D}(order_high::Int, order_low::Int) where {D}
        @assert isodd(order_high) && isodd(order_low)
        @assert order_high > order_low ≥ 1
        return new{D}(order_high, order_low)
    end

    function GrundmannMoeller{D}(order_high::Int) where {D}
        @assert isodd(order_high)
        @assert order_high ≥ 3
        return new{D}(order_high, order_high - 2)
    end
end
