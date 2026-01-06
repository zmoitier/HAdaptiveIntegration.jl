"""
    @kwdef struct TabulatedEmbeddedCubature{DOM<:AbstractDomain} <: AbstractRule{DOM}

An embedded cubature rule consisting of a high order cubature rule and a low order cubature
rule. Note that the low order cubature uses `nodes[1:L]` as its nodes where `L` is the
length of the `weights_low`. The cubature nodes and weights are assumed to be for the
reference domain (use the [`reference_domain`](@ref) function to get the reference domain).

## Fields:
- `description::String`: description of the embedded cubature.
- `reference::String`: where the values are found.
- `order_high::Int`: order of the high order cubature.
- `order_low::Int`: order of the low order cubature.
- `precision::Int`: number of significant digits on the node and weight values,
  `10^-precision` give the relative precision of the values.
- `nodes::Vector{Vector{String}}`: the cubature nodes.
- `weights_high::Vector{String}`: the cubature weights for the high order cubature.
- `weights_low::Vector{String}`: the cubature weights for the low order cubature.

## Invariants (check at construction):
- `length(nodes) == length(weights_high)`
- `length(weights_high) ≥ length(weights_low)`
- `order_high ≥ order_low` ≥ 0
- `precision ≥ 0`
"""
struct TabulatedEmbeddedCubature{DOM<:AbstractDomain} <: AbstractRule{DOM}
    description::String
    reference::String
    order_high::Int
    order_low::Int
    precision::Int
    nodes::Vector{Vector{String}}
    weights_high::Vector{String}
    weights_low::Vector{String}

    function TabulatedEmbeddedCubature{DOM}(;
        description::String,
        reference::String,
        order_high::Int,
        order_low::Int,
        precision::Int,
        nodes::Vector{Vector{String}},
        weights_high::Vector{String},
        weights_low::Vector{String},
    ) where {DOM<:AbstractDomain}
        D = dimension(DOM)
        @assert all(n -> length(n) == D, nodes) "Each node must have length equal to the \
        dimension D"
        @assert length(nodes) == length(weights_high) "The number of nodes must match the \
        number of high-order weights"
        @assert length(weights_high) ≥ length(weights_low) "The length of high order \
        weights must be greater or equal to the length of low-order weights"
        @assert order_high ≥ order_low ≥ 0 "order_high must be greater than or equal to \
        order_low and orders must be non-negative"
        @assert precision ≥ 0
        return new{DOM}(
            description,
            reference,
            order_high,
            order_low,
            precision,
            nodes,
            weights_high,
            weights_low,
        )
    end
end

"""
    orders(rule::AR) where {AR<:AbstractRule}

Return the high and low order of the embedded cubature `rule`.
"""
function orders(tec::TabulatedEmbeddedCubature)
    return tec.order_high, tec.order_low
end

function embedded_cubature(
    tec::TabulatedEmbeddedCubature{DOM}, (::Type{T})=float(Int)
) where {DOM<:AbstractDomain,T}
    if eps(T) < 10.0^(-tec.precision)
        @warn "The embedded cubature `$(tec.description)` has fewer significant digits \
        than type $T, which may lead to inaccurate computations."
    end
    D = dimension(DOM)
    return EmbeddedCubature(
        [SVector{D,T}(parse.(T, x)) for x in tec.nodes],
        parse.(T, tec.weights_high),
        parse.(T, tec.weights_low),
    )
end
