"""
    struct EmbeddedCubature{D,T}

An embedded cubature rule consisting of a high order cubature rule with `H` nodes and a low
order cubature rule with `L` nodes. Note that the low order cubature uses `nodes[1:L]` as
its nodes. The cubature nodes and weights are assume to be for the reference domain.

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
        @assert length(nodes) == length(weights_high)
        @assert length(weights_high) ≥ length(weights_low)
        return new{D,T}(nodes, weights_high, weights_low)
    end
end

"""
    embedded_cubature(tec::TabulatedEmbeddedCubature)
    embedded_cubature(T::DataType, tec::TabulatedEmbeddedCubature)

    embedded_cubature(gm::GrundmannMoeller)
    embedded_cubature(T::DataType, gm::GrundmannMoeller)

    embedded_cubature(nodes, weights_high, weights_low)
    embedded_cubature(T::DataType, nodes, weights_high, weights_low)

Construct an embedded cubature from a [`TabulatedEmbeddedCubature`](@ref), a
[`GrundmannMoeller`](@ref), cubature, or from a vector of nodes and two vectors of weights
for the high and low order cubature, with element type `T`.
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

function embedded_cubature(T::DataType, tec::TabulatedEmbeddedCubature)
    if 10 * eps(T) < 10.0^(-tec.nb_significant_digits)
        @warn "the embedded cubature `$(tec.description)` has less precision than type $T."
    end

    return EmbeddedCubature(
        [SVector(parse.(T, x)) for x in tec.nodes],
        parse.(T, tec.weights_high),
        parse.(T, tec.weights_low),
    )
end
embedded_cubature(tec::TabulatedEmbeddedCubature) = embedded_cubature(float(Int), tec)

function embedded_cubature(T::DataType, gm::GrundmannMoeller)
    vol = 1 / T(factorial(gm.dim)) # volume of the reference simplex

    # high order cubature
    Tn = grundmann_moeller(T, Val(gm.dim), gm.deg)

    H = length(Tn.points)
    nodes_high = [SVector{gm.dim}(Tn.points[H - i + 1][2:end]) for i in 1:H]
    weights_high = [Tn.weights[H - i + 1] * vol for i in 1:H]

    # low order cubature
    Tn_low = grundmann_moeller(T, Val(gm.dim), gm.deg - 2)
    L = length(Tn_low.points)
    weights_low = [Tn_low.weights[L - i + 1] * vol for i in 1:L]

    return EmbeddedCubature(nodes_high, weights_high, weights_low)
end
embedded_cubature(gm::GrundmannMoeller) = embedded_cubature(float(Int), gm)

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
    fct, domain::AbstractDomain{D}, norm=x -> LinearAlgebra.norm(x, Inf)
) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)
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
