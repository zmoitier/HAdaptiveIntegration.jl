module IncreasePrecisionExt

import HAdaptiveIntegration: increase_precision

using Base.Iterators: countfrom, partition
using ForwardDiff.DiffResults: JacobianResult, MutableDiffResult, jacobian, value
using ForwardDiff: jacobian!
using HAdaptiveIntegration:
    Orthotope,
    Segment,
    Simplex,
    TabulatedEmbeddedCubature,
    dimension,
    integral_monomials_exact,
    orders
using LinearAlgebra: norm
using Printf: @sprintf
using StaticArrays: SVector

function __init__()
    @info "Loading IncreasePrecisionExt.jl"
end

mutable struct NewtonState{T}
    iter::Int
    x::Vector{T}
    result::MutableDiffResult{1,Vector{T},Tuple{Matrix{T}}}
    Δ::Vector{T}
    x_abs_err::T
    f_abs_err::T

    function NewtonState(x₀::Vector{T}, f::Function) where {T}
        result = JacobianResult(f(x₀), x₀)
        jacobian!(result, f, x₀)
        Δ = jacobian(result) \ value(result)
        x₁ = x₀ - Δ
        x_abs_err = norm(Δ, Inf)
        f_abs_err = norm(value(result), Inf)
        return new{T}(1, x₁, result, Δ, x_abs_err, f_abs_err)
    end
end

function iter!(ns::NewtonState{T}, f::Function) where {T}
    jacobian!(ns.result, f, ns.x)
    ns.Δ .= jacobian(ns.result) \ value(ns.result)
    ns.x .-= ns.Δ
    ns.x_abs_err = norm(ns.Δ, Inf)
    ns.f_abs_err = norm(value(ns.result), Inf)
    ns.iter += 1
    return ns
end

function newton(f, x₀, x_atol, f_atol, maxiter)
    ns = NewtonState(x₀, f)

    @info @sprintf(
        "iter %3d: |xₙ - xₙ₋₁| = %.2e, |f(xₙ)| = %.2e", ns.iter, ns.x_abs_err, ns.f_abs_err,
    )

    while (ns.x_abs_err > x_atol) && (ns.f_abs_err > f_atol) && (ns.iter < maxiter)
        iter!(ns, f)
        @info @sprintf(
            "iter %3d: |xₙ - xₙ₋₁| = %.2e, |f(xₙ)| = %.2e",
            ns.iter,
            ns.x_abs_err,
            ns.f_abs_err,
        )
    end

    @info """Newton convergence report:

          iteration = $(ns.iter) / $maxiter
        |xₙ - xₙ₋₁| = $(@sprintf("%.2e", ns.x_abs_err)) $(ns.x_abs_err > x_atol ? "≰" : "≤") $(@sprintf("%.2e", x_atol))
            |f(xₙ)| = $(@sprintf("%.2e", ns.f_abs_err)) $(ns.f_abs_err > f_atol ? "≰" : "≤") $(@sprintf("%.2e", f_atol))
    """

    if ns.iter ≥ maxiter
        @warn "maximum number of iterations reached, try increasing the keyword argument `maxiter=$maxiter`."
    end

    return ns.x, ns.x_abs_err
end

function increase_precision(
    tec::TabulatedEmbeddedCubature{DOM},
    ::Type{T};
    x_atol=10 * eps(T),
    f_atol=10 * eps(T),
    maxiter::Int=16,
) where {DOM,T}
    if eps(T) ≥ 10.0^(-tec.precision)
        @info "Tabulated precision $(tec.precision) is sufficient for type $(T). No need to increase precision."
        return tec
    end

    @info "Increasing to target precision $(-floor(Int, log10(x_atol))) from tabulated precision $(tec.precision). This may take some time."

    D = dimension(DOM)
    U, range_nodes, range_wh, range_wl = pack(tec, T, D)
    order_high, order_low = orders(tec)
    exponent2values, range_low, range_high = _flatten_with_type(
        integral_monomials_exact(DOM, order_high), T, order_high, order_low
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

    U, δ = newton(G, U, x_atol, f_atol, maxiter)

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
