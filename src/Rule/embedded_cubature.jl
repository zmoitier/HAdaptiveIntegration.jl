"""
    struct EmbeddedCubature{D,T}

An embedded cubature rule consisting of a high order cubature rule with `H` nodes and a low
order cubature rule with `L` nodes. Note that the low order cubature uses `nodes[1:L]` as
its nodes. The cubature nodes and weights are assumed to be for the reference domain.

## Fields:
- `nodes::Vector{SVector{D,T}}`: the cubature nodes.
- `weights_high::Vector{T}`: the cubature weights for the high order cubature.
- `weights_low::Vector{T}`: the cubature weights for the low order cubature.

## Invariants (check at construction):
- `length(nodes) == length(weights_high)`
- `length(weights_high) ≥ length(weights_low)`
"""
struct EmbeddedCubature{D,T}
    nodes::Vector{SVector{D,T}}
    weights_high::Vector{T}
    weights_low::Vector{T}

    function EmbeddedCubature(
        nodes::Vector{SVector{D,T}}, weights_high::Vector{T}, weights_low::Vector{T}
    ) where {D,T}
        @assert length(nodes) == length(weights_high) "The number of nodes must match the number of high-order weights."
        @assert length(weights_high) ≥ length(weights_low) "weights_high must have a length greater than or equal to weights_low."
        return new{D,T}(nodes, weights_high, weights_low)
    end
end

"""
    embedded_cubature(ar::AbstractRule)
    embedded_cubature(T::DataType, ar::AbstractRule)

    embedded_cubature(nodes, weights_high, weights_low)
    embedded_cubature(T::DataType, nodes, weights_high, weights_low)

Construct the embedded cubature with element type `T` from an `AbstractRule`
([`TabulatedEmbeddedCubature`](@ref) or [`GrundmannMoeller`](@ref)), or from a vector of
nodes and two vectors of weights for the high and low order cubature, with element type `T`.
"""
function embedded_cubature(T::DataType, nodes, weights_high, weights_low)
    @assert allequal(length, nodes) "all nodes should have the same length."
    D = length(first(nodes))

    return EmbeddedCubature(SVector{D,T}.(nodes), weights_high, weights_low)
end
function embedded_cubature(nodes, weights_high, weights_low)
    T_nodes = promote_to_float(nodes...)
    T_weights = promote_to_float(weights_high, weights_low)
    T = promote_type(T_nodes, T_weights)
    return embedded_cubature(float(T), nodes, weights_high, weights_low)
end

"""
    (ec::EmbeddedCubature{D,T})(
        fct, domain::Domain{D,T}, norm=x -> LinearAlgebra.norm(x, Inf)
    ) where {D,T}

Return `I_high` and `norm(I_high - I_low)` where `I_high` and `I_low` are the result of the
high order cubature and the low order cubature on `domain`. The function `fct` must take a
`SVector{D,T}` to a return type `K`, and `K` must support the multiplication by a scalar and
the addition. Note that there is no check, beyond compatibility of dimension and type, that
the embedded cubature is for the right domain.
"""
function (ec::EmbeddedCubature{D,T})(
    fct, domain::AbstractDomain{D,T}, norm=x -> norm(x, Inf)
) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)
    μ = abs_det_jac(domain)
    Φ = map_from_reference(domain)

    v = fct(Φ(ec.nodes[1]))
    Iₗ = ec.weights_low[1] * v
    Iₕ = ec.weights_high[1] * v
    for i in 2:L
        v = fct(Φ(ec.nodes[i]))
        Iₗ += ec.weights_low[i] * v
        Iₕ += ec.weights_high[i] * v
    end

    for i in (L + 1):H
        Iₕ += ec.weights_high[i] * fct(Φ(ec.nodes[i]))
    end

    return μ * Iₕ, μ * norm(Iₕ - Iₗ)
end
