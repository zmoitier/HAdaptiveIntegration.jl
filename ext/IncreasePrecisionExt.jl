module IncreasePrecisionExt

import HAdaptiveIntegration: increase_precision

using Base.Iterators: countfrom, partition
using ForwardDiff: jacobian
using HAdaptiveIntegration:
    Orthotope, Segment, Simplex, TabulatedEmbeddedCubature, dimension, orders
using LinearAlgebra: norm
using StaticArrays: SVector

function __init__()
    @info "Loading IncreasePrecisionExt.jl"
end

function increase_precision(
    tec::TabulatedEmbeddedCubature{DOM}, ::Type{T}; atol::T=10 * eps(T)
) where {DOM,T}
    if eps(T) ≥ 10.0^(-tec.precision)
        @info "Tabulated precision $(tec.precision) is sufficient for type $T. No need to increase precision."
        return tec
    end
    @info "Need to increase precision, it may take some time."

    D = dimension(DOM)
    U, range_nodes, range_wh, range_wl = pack(tec, T, D)
    order_high, order_low = orders(tec)

    if DOM <: Segment
        exp2int = exp2int = [[(i,) => 1//(i + 1)] for i in 0:order_high]
    elseif DOM <: Simplex
        exp2int = integral_monomial_simplex(D, order_high)
    elseif DOM <: Orthotope
        exp2int = integral_monomial_orthotope(D, order_high)
    else
        error("increasing_precision is not implemented for domain type: $DOM")
    end

    L = length(range_wl)
    function G(U::Vector{S}) where {S}
        nodes = @view U[range_nodes]
        wh = @view U[range_wh]
        wl = @view U[range_wl]

        V = Vector{S}(undef, order_high + order_low + 2)
        for (k, e2i) in enumerate(exp2int[1:(order_low + 1)])
            for (e, v) in e2i
                vh, vl = zero(S), zero(S)
                for (i, node) in enumerate(partition(nodes, D))
                    p = prod(node .^ e)
                    vh += wh[i] * p
                    if i ≤ L
                        vl += wl[i] * p
                    end
                end
                V[k] = vh - v
                V[order_high + 1 + k] = vl - v
            end
        end
        for (k, e2i) in enumerate(exp2int[(order_low + 2):(order_high + 1)])
            for (e, v) in e2i
                vh = zero(v)
                for (i, node) in enumerate(partition(nodes, D))
                    p = prod(node .^ e)
                    vh += wh[i] * p
                end
                V[order_low + 1 + k] = vh - v
            end
        end
        return V
    end

    U, δ = newton!(U, G, x -> jacobian(G, x); atol=atol)

    nodes, weights_high, weights_low = unpack(U, range_nodes, range_wh, range_wl, D)
    return TabulatedEmbeddedCubature{DOM}(;
        description=tec.description * " (increased precision)",
        reference=tec.reference,
        order_high=tec.order_high,
        order_low=tec.order_low,
        precision=-floor(Int, log10(δ)),
        nodes=[Vector(string.(node)) for node in nodes],
        weights_high=string.(weights_high),
        weights_low=string.(weights_low),
    )
end

function newton!(U, G, DG; norm=x -> norm(x, Inf), atol=10 * eps(T), maxiter::Int=16)
    iter = 0
    R = G(U)
    Δ = DG(U) \ R
    U = U - Δ

    iter = 1
    while (norm(Δ) > atol) && (norm(R) > atol) && (iter < maxiter)
        R = G(U)
        Δ = DG(U) \ R
        U = U - Δ
        iter += 1
    end

    @info """Newton iteration : $iter

        |Uₙ - Uₙ₋₁| = $(norm(Δ))
            |F(Uₙ)| = $(norm(R))
    """

    if iter ≥ maxiter
        @warn "maximum number of iterations reached, try increasing the keyword argument `maxiter=$maxiter`."
    end

    return U, norm(Δ)
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

function integral_monomial_simplex(dim::Int, deg_tot_max::Int)
    @assert (dim > 0) && (deg_tot_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    exp2int = [[(i,) => 1//prod((i + 1):(i + dim))] for i in 0:deg_tot_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, e2i) in zip(countfrom(0), exp2int)
            for (e, v) in e2i
                push!(new[k + 1], (0, e...) => v)
                t = 1
                for n in 1:(deg_tot_max - k)
                    t *= n//(n + k + dim)
                    push!(new[k + 1 + n], (n, e...) => t * v)
                end
            end
        end
        exp2int = new
    end

    return exp2int
end

function integral_monomial_orthotope(dim::Int, deg_tot_max::Int)
    @assert (dim > 0) && (deg_tot_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    exp2int = [[(i,) => 1//(i + 1)] for i in 0:deg_tot_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, e2i) in zip(countfrom(0), exp2int)
            for (e, v) in e2i
                for n in 0:(deg_tot_max - k)
                    push!(new[k + 1 + n], (n, e...) => v//(n + 1))
                end
            end
        end
        exp2int = new
    end

    return exp2int
end

end
