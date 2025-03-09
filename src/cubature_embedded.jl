"""
    struct TabulatedEmbeddedCubature

An embedded cubature rule consisting of a high order cubature rule and a low order cubature rule.
Note that the low order cubature uses `nodes[1:L]` as its nodes where `L` is the length of the `weights_low`.
The cubature nodes and weights are assume to be for the reference domain.

## Fields:
- `name::String`: name of the embedded cubature.
- `domain::String`: domain of the cubature, the possible values are "reference segment", "reference rectangle", "reference cuboid", "reference triangle", and "reference tetrahedron".
- `reference::String`: where the values are found.
- `nb_significant_digits::Int`: number of significant digits on the node and weight values, `10^-nb_significant_digits` give the relative precision of the values.
- `nodes::Vector{Vector{String}}`: the cubature nodes.
- `weights_high::Vector{String}`: the cubature weights for the high order cubature.
- `order_high::Int`: order of the high order cubature.
- `weights_low::Vector{String}`: the cubature weights for the low order cubature.
- `order_low::Int`: order of the low order cubature.
"""
@kwdef struct TabulatedEmbeddedCubature
    name::String
    domain::String
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
Note that the low order cubature uses `nodes[1:L]` as its nodes.
The cubature nodes and weights are assume to be for the reference domain.

## Fields:
- `nodes::SVector{H,SVector{D,T}}`: the cubature nodes.
- `weights_high::SVector{H,T}`: the cubature weights for the high order cubature.
- `weights_low::SVector{L,T}`: the cubature weights for the low order cubature.
"""
struct EmbeddedCubature{H,L,D,T}
    nodes::SVector{H,SVector{D,T}}
    weights_high::SVector{H,T}
    weights_low::SVector{L,T}
end

"""
    embedded_cubature(
        nodes::Vector{VT}, weights_high::VT, weights_low::VT
    ) where {VT<:AbstractVector{<:Real}}

Return an embedded cubature form a vector of nodes and two vector of weights for the high order and low order cubature.
"""
function embedded_cubature(
    nodes::Vector{VT}, weights_high::VT, weights_low::VT
) where {VT<:AbstractVector{<:Real}}
    @assert allequal(length, nodes) "all nodes should have the same length."
    D = length(first(nodes))

    @assert length(nodes) == length(weights_high) "nodes and weights_high should have the same length."
    H = length(weights_high)

    @assert length(weights_high) ≥ length(weights_low) "weights_high length should be greater than weights_low length."
    L = length(weights_low)

    return EmbeddedCubature(
        SVector{H}(SVector{D}.(nodes)), SVector{H}(weights_high), SVector{L}(weights_low)
    )
end

"""
    embedded_cubature(tec::TabulatedEmbeddedCubature, T::DataType=Float64)

Return the embedded cubature with value type `T` from a `TabulatedEmbeddedCubature`.
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

"""
   struct GrundmannMoeller

Cubature rule for a `dim`-simplex of degree `deg`.
"""
struct GrundmannMoeller
    dim::Int
    deg::Int
end

"""
    embedded_cubature(gm::GrundmannMoeller, T::DataType=Float64)

Return the embedded cubature of Grundmann and Möller.
"""
function embedded_cubature(gm::GrundmannMoeller, T::DataType=Float64)
    dimension, degree = gm.dim, gm.deg
    vol = 1 / T(factorial(dimension)) # volume of the reference simplex

    # high order cubature
    Tn = grundmann_moeller(T, Val(dimension), degree)
    H = length(Tn.points)
    nodes_high = SVector(ntuple(i -> SVector{dimension}(Tn.points[H - i + 1][2:end]), H))
    weights_high = SVector(ntuple(i -> Tn.weights[H - i + 1] * vol, H))

    # low order cubature
    Tn_low = grundmann_moeller(T, Val(dimension), degree - 2)
    L = length(Tn_low.points)
    weights_low = SVector(ntuple(i -> Tn_low.weights[L - i + 1] * vol, L))

    return EmbeddedCubature(nodes_high, weights_high, weights_low)
end

"""
    (ec::EmbeddedCubature{H,L,D,T})(
        fct, domain::Domain{D,T}, norm=x -> LinearAlgebra.norm(x, Inf)
    ) where {H,L,D,T}

Return `I_high` and `norm(I_high - I_low)` where `I_high` and `I_low` are the result of the high order cubature and the low order cubature on `domain`.
The function `fct` must take a `SVector{D,T}` to a return type `K`, and `K` must support the multiplication by a scalar and the addition.
Note that there is no check, beyond compatibility of dimension and type, that the embedded cubature is for the right domain.
"""
function (ec::EmbeddedCubature{H,L,D,T})(
    fct, domain::Domain{D,T}, norm=x -> LinearAlgebra.norm(x, Inf)
) where {H,L,D,T}
    μ = abs_det_jac(domain)
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

    return μ * I_high, μ * norm(I_high - I_low)
end
