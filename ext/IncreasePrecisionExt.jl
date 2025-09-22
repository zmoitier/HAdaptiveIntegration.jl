module IncreasePrecisionExt

import HAdaptiveIntegration: increase_precision

using Base.Iterators: countfrom, partition
using ForwardDiff: jacobian
using HAdaptiveIntegration:
    Orthotope,
    Segment,
    Simplex,
    TabulatedEmbeddedCubature,
    dimension,
    integral_monomials,
    orders
using LinearAlgebra: norm
using Printf: @sprintf
using StaticArrays: SVector

function __init__()
    @info "Loading IncreasePrecisionExt.jl"
end

mutable struct NewtonState{T}
    iter::Int
    U::Vector{T}
    R::Vector{T}
    r::T
    J::Matrix{T}
    Δ::Vector{T}
    δ::T
end

function increase_precision(
    tec::TabulatedEmbeddedCubature{DOM}, ::Type{T}; atol::T=10 * eps(T), maxiter::Int=16
) where {DOM,T}
    if eps(T) ≥ 10.0^(-tec.precision)
        @info LazyString(
            "Tabulated precision ",
            tec.precision,
            " is sufficient for type ",
            T,
            ". No need to increase precision.",
        )
        return tec
    end

    @info LazyString(
        "Increasing to target precision ",
        -floor(Int, log10(atol)),
        " from tabulated precision ",
        tec.precision,
        ". This may take some time.",
    )

    D = dimension(DOM)
    U, range_nodes, range_wh, range_wl = pack(tec, T, D)
    order_high, order_low = orders(tec)
    exponent2values, range_low, range_high = _flatten_with_type(
        integral_monomials(DOM, order_high), T, order_high, order_low
    )

    L = length(range_wl)
    function G(U::Vector{S}) where {S}
        nodes = @view U[range_nodes]
        wh = @view U[range_wh]
        wl = @view U[range_wl]

        V = Vector{S}(undef, range_high.stop + range_low.stop)
        for (k, (α, vₑₓ)) in enumerate(exponent2values[range_low])
            vh, vl = zero(S), zero(S)
            for (i, node) in enumerate(partition(nodes, D))
                r = prod(node .^ α)
                vh += wh[i] * r
                if i ≤ L
                    vl += wl[i] * r
                end
            end
            V[k] = vh - vₑₓ
            V[range_high.stop + k] = vl - vₑₓ
        end
        for (k, (α, vₑₓ)) in enumerate(exponent2values[range_high])
            vh = zero(S)
            for (i, node) in enumerate(partition(nodes, D))
                vh += wh[i] * prod(node .^ α)
            end
            V[range_low.stop + k] = vh - vₑₓ
        end

        return V
    end

    U, δ = newton(G, x -> jacobian(G, x), U, atol, maxiter)

    nodes, weights_high, weights_low = unpack(U, range_nodes, range_wh, range_wl, D)
    return TabulatedEmbeddedCubature{DOM}(;
        description=(tec.description * " (increased precision)"),
        reference=tec.reference,
        order_high=tec.order_high,
        order_low=tec.order_low,
        precision=(-floor(Int, log10(δ))),
        nodes=[Vector(string.(node)) for node in nodes],
        weights_high=string.(weights_high),
        weights_low=string.(weights_low),
    )
end

function newton(G, DG, U, atol, maxiter)
    V = copy(U)

    iter = 0
    R = G(V)
    r = norm(R, Inf)
    Δ = DG(V) \ R
    δ = norm(Δ, Inf)
    V = V - Δ

    iter = 1
    while (δ > atol) && (r > atol) && (iter < maxiter)
        R = G(V)
        r = norm(R, Inf)
        Δ = DG(V) \ R
        @show δ = norm(Δ, Inf)
        V = V - Δ
        iter += 1
    end

    @info """Newton convergence report:

    |Uₙ - Uₙ₋₁| = $(@sprintf("%.2e", δ)) $(δ > atol ? "≰" : "≤") $(@sprintf("%.2e", atol))
        |F(Uₙ)| = $(@sprintf("%.2e", r)) $(r > atol ? "≰" : "≤") $(@sprintf("%.2e", atol))
     #iteration = $iter $(iter < maxiter ? "<" : "≮") $maxiter
    """

    if iter ≥ maxiter
        @warn "maximum number of iterations reached, try increasing the keyword argument `maxiter=$maxiter`."
    end

    return V, δ
end

# [ nodes[1]..., ..., nodes[end]..., weights_high..., weights_low... ]
function pack(tec::TabulatedEmbeddedCubature, ::Type{T}, D::Int) where {T}
    U = Vector{T}()
    for node in tec.nodes
        append!(U, parse.(T, node))
    end
    append!(U, parse.(T, tec.weights_high))
    append!(U, parse.(T, tec.weights_low))

    H, L = length(tec.weights_high), length(tec.weights_low)
    return U, 1:(D * H), (H * D + 1):(H * D + H), (H * D + H + 1):(H * D + H + L)
end

function unpack(
    U::Vector{T},
    range_pts::UnitRange{Int},
    range_wh::UnitRange{Int},
    range_wl::UnitRange{Int},
    D::Int,
) where {T}
    nodes = Vector{SVector{D,T}}()
    for node in partition(U[range_pts], D)
        push!(nodes, SVector{D}(node))
    end

    weights_high = U[range_wh]
    weights_low = U[range_wl]

    return (nodes, weights_high, weights_low)
end

function _flatten_with_type(
    exponent2values::Vector{Vector{Pair{NTuple{D,Int},S}}},
    ::Type{T},
    order_high::Int,
    order_low::Int,
) where {D,S,T}
    nb_lo, nb_hi = binomial(order_low + D, D), binomial(order_high + D, D)

    result = Vector{Pair{NTuple{D,Int},T}}()
    for α2v in exponent2values
        append!(result, α2v)
    end

    return result, 1:nb_lo, (nb_lo + 1):nb_hi
end

end
