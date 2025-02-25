"""
    struct EmbeddedCubatureRaw

An embedded cubature rule consisting of a high order cubature rule nodes and a low order cubature rule.
The cubature nodes and weights are assume to be for the reference simplex or orthotope.
Note that the low order cubature uses `nodes[1:L]` as its nodes.

## Fields:
- `name::String`: name of the embedded cubature;
- `reference::String`: where the values are found;
- `nb_significant_digits::Int`: number of significant digits on the node and weight values, `10^-nb_significant_digits` give the relative precision of the values;
- `nodes::SVector{H,SVector{D,T}}`: the cubature nodes;
- `weights_high::SVector{H,T}`: the cubature weights for the high order cubature;
- `order_high::Int`: order of the high order cubature;
- `weights_low::SVector{L,T}`: the cubature weights for the low order cubature;
- `order_low::Int`: order of the low order cubature.
"""
@kwdef struct EmbeddedCubatureRaw
    name::String
    reference::String
    nb_significant_digits::Int
    nodes::Vector{Vector{String}}
    weights_high::Vector{String}
    order_high::Int
    weights_low::Vector{String}
    order_low::Int
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

Construct an embedded cubature form a vector of nodes and two vector of weights high and low.
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

function embedded_cubature_from_raw(ecr::EmbeddedCubatureRaw, T::DataType=Float64)
    @assert eps(T) > 10.0^(-ecr.nb_significant_digits) "the embedded cubature `$(ecr.name)` has less precision than type $T."

    return embedded_cubature(
        [parse.(T, x) for x in ecr.nodes],
        parse.(T, ecr.weights_high),
        parse.(T, ecr.weights_low),
    )
end

function (ec::EmbeddedCubature{H,L,D,T})(
    fct::Function, norm=LinearAlgebra.norm
) where {H,L,D,T<:Real}
    v = fct(ec.nodes[1])
    I_low = v * ec.weights_low[1]
    I_high = v * ec.weights_high[1]
    for i in 2:L
        v = fct(ec.nodes[i])
        I_low += v * ec.weights_low[i]
        I_high += v * ec.weights_high[i]
    end

    for i in (L + 1):H
        I_high += fct(ec.nodes[i]) * ec.weights_high[i]
    end

    # TODO: check if it is better to use relative error.
    return I_high, norm(I_high - I_low, Inf)
end

function (ec::EmbeddedCubature{H,L,D,T})(
    fct::Function, domain::Domain{D,T}, norm=LinearAlgebra.norm
) where {H,L,D,T<:Real}
    μ = abs_jacobian_determinant(domain)
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

    # TODO: check if it is better to use relative error.
    return μ * I_high, μ * norm(I_high - I_low, Inf)
end
