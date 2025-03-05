"""
    struct TabulatedEmbeddedCubature

An embedded cubature rule consisting of a high order cubature rule nodes and a low order cubature rule.
The cubature nodes and weights are assume to be for the reference simplex or orthotope.
Note that the low order cubature uses `nodes[1:L]` as its nodes.

## Fields:
- `name::String`: name of the embedded cubature;
- `domain::Domain`: domain of the cubature.
- `reference::String`: where the values are found;
- `nb_significant_digits::Int`: number of significant digits on the node and weight values, `10^-nb_significant_digits` give the relative precision of the values;
- `nodes::Vector{Vector{String}}`: the cubature nodes;
- `weights_high::Vector{String}`: the cubature weights for the high order cubature;
- `order_high::Int`: order of the high order cubature;
- `weights_low::Vector{String}`: the cubature weights for the low order cubature;
- `order_low::Int`: order of the low order cubature.
"""
struct TabulatedEmbeddedCubature
    name::String
    domain::Domain
    reference::String
    nb_significant_digits::Int
    nodes::Vector{Vector{String}}
    weights_high::Vector{String}
    order_high::Int
    weights_low::Vector{String}
    order_low::Int
end

"""
    tabulated_embedded_cubature(;
        name::String,
        domain::Domain,
        reference::String,
        nb_significant_digits::Int,
        nodes::Vector{Vector{String}},
        weights_high::Vector{String},
        order_high::Int,
        weights_low::Vector{String},
        order_low::Int,
    )

An embedded cubature rule consisting of a high order cubature rule nodes and a low order cubature rule.
The cubature nodes and weights are assume to be for the reference simplex or orthotope.
Note that the low order cubature uses `nodes[1:L]` as its nodes.

## Fields:
- `name::String`: name of the embedded cubature;
- `domain::Domain`: domain of the cubature.
- `reference::String`: where the values are found;
- `nb_significant_digits::Int`: number of significant digits on the node and weight values, `10^-nb_significant_digits` give the relative precision of the values;
- `nodes::Vector{Vector{String}}`: the cubature nodes;
- `weights_high::Vector{String}`: the cubature weights for the high order cubature;
- `order_high::Int`: order of the high order cubature;
- `weights_low::Vector{String}`: the cubature weights for the low order cubature;
- `order_low::Int`: order of the low order cubature.
"""
function tabulated_embedded_cubature(;
    name::String,
    domain::Domain{D,T},
    reference::String,
    nb_significant_digits::Int,
    nodes::Vector{Vector{String}},
    weights_high::Vector{String},
    order_high::Int,
    weights_low::Vector{String},
    order_low::Int,
) where {D,T<:Real}
    @assert all(length(node) == D for node in nodes) "all nodes must have length `D` where `D` is the domain dimension."

    return TabulatedEmbeddedCubature(
        name,
        domain,
        reference,
        nb_significant_digits,
        nodes,
        weights_high,
        order_high,
        weights_low,
        order_low,
    )
end

"""
    struct EmbeddedCubature{H,L,D,T}

An embedded cubature rule consisting of a high order cubature rule with `H` nodes and a low order cubature rule with `L` nodes.
The cubature nodes and weights are assume to be for the reference simplex or orthotope.
Note that the low order cubature uses `nodes[1:L]` as its nodes.

## Fields:
- `nodes::SVector{H,SVector{D,T}}`: the cubature nodes;
- `weights_high::SVector{H,T}`: the cubature weights for the high order cubature;
- `weights_low::SVector{L,T}`: the cubature weights for the low order cubature.
"""
struct EmbeddedCubature{H,L,D,T}
    nodes::SVector{H,SVector{D,T}}
    weights_high::SVector{H,T}
    weights_low::SVector{L,T}
end

"""
    embedded_cubature(nodes, weights_high, weights_low)

Return an embedded cubature form a vector of nodes and two vector of weights for the high order and low order cubature.
"""
function embedded_cubature(
    nodes::Vector{Vector{T}}, weights_high::Vector{T}, weights_low::Vector{T}
) where {T<:Real}
    @assert allequal(length, nodes) "all nodes should have the same length."
    D = length(first(nodes))

    @assert length(nodes) == length(weights_high) "nodes and weights_high should have the same length."
    H = length(weights_high)

    @assert length(weights_high) > length(weights_low) "weights_high length should be greater than weights_low length."
    L = length(weights_low)

    return EmbeddedCubature(
        SVector{H,SVector{D,T}}(SVector{D,T}.(nodes)),
        SVector{H,T}(weights_high),
        SVector{L,T}(weights_low),
    )
end

"""
    embedded_cubature(tec::TabulatedEmbeddedCubature, T=Float64)

Return the `EmbeddedCubature` with type `T` from an `TabulatedEmbeddedCubature`.
"""
function embedded_cubature(tec::TabulatedEmbeddedCubature, T::DataType=Float64)
    if 10 * eps(T) < 10.0^(-tec.nb_significant_digits)
        @warn "the embedded cubature `$(tec.name)` has less precision than type $T."
    end

    return embedded_cubature(
        [parse.(T, x) for x in tec.nodes],
        parse.(T, tec.weights_high),
        parse.(T, tec.weights_low),
    )
end

function (ec::EmbeddedCubature{H,L,D,T})(
    fct::Function, domain::Domain{D,T}, norm=LinearAlgebra.norm
) where {H,L,D,T<:Real}
    μ = abs_det_jacobian(domain)
    Φ = map_from_reference(domain)

    v = fct(Φ(ec.nodes[1]))
    I_low = v * ec.weights_low[1]
    I_high = v * ec.weights_high[1]
    for i in 2:L
        v = fct(Φ(ec.nodes[i]))
        I_low += v * ec.weights_low[i]
        I_high += v * ec.weights_high[i]
    end

    for i in (L + 1):H
        I_high += fct(Φ(ec.nodes[i])) * ec.weights_high[i]
    end

    return μ * I_high, μ * norm(I_high - I_low, Inf)
end
