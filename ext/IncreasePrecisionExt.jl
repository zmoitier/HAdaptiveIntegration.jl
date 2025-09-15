module IncreasePrecisionExt

import HAdaptiveIntegration: increase_precision

using Base.Iterators: countfrom, partition
using HAdaptiveIntegration:
    EmbeddedCubature,
    Orthotope,
    Simplex,
    TabulatedEmbeddedCubature,
    dimension,
    embedded_cubature,
    orders
using Optim: LBFGS, Options, minimizer, optimize
using StaticArrays: SVector

function __init__()
    @info "Loading IncreasePrecisionExt.jl"
end

function increase_precision(
    tec::TabulatedEmbeddedCubature{DOM}, ::Type{T}, options::Options=Options()
) where {DOM,T}
    εT = eps(T)

    if εT ≥ 10.0^(-tec.precision)
        @info "Tabulated precision $(tec.precision) is sufficient for type $T. No need to increase precision."
        return embedded_cubature(tec, T)
    end
    @info "Need to increase precision, it may take some time."

    D = dimension(DOM)
    U, range_nodes, range_wh, range_wl = pack(tec, T, D)
    order_high, order_low = orders(tec)

    if DOM <: Simplex
        exp2int = integral_monomial_simplex(D, order_high)
    elseif DOM <: Orthotope
        exp2int = integral_monomial_orthotope(D, order_high)
    else
        error("increasing_precision is not implemented for domain type: $DOM")
    end

    L = length(range_wl)
    function F(U)
        r = zero(eltype(U))

        nodes = @view U[range_nodes]
        wh = @view U[range_wh]
        wl = @view U[range_wl]

        for e2i in exp2int[1:order_low]
            for (e, v) in e2i
                vh, vl = zero(v), zero(v)
                for (i, node) in enumerate(partition(nodes, D))
                    p = prod(node .^ e)
                    vh += wh[i] * p
                    if i ≤ L
                        vl += wl[i] * p
                    end
                end
                r += (vh - v)^2 + (vl - v)^2
            end
        end
        for e2i in exp2int[(order_low + 1):order_high]
            for (e, v) in e2i
                vh = zero(v)
                for (i, node) in enumerate(partition(nodes, D))
                    p = prod(node .^ e)
                    vh += wh[i] * p
                end
                r += (vh - v)^2
            end
        end

        return r
    end

    @info "F(U₀) = $(F(U))"

    result = optimize(F, U, LBFGS(), options; autodiff=:forward)
    @info result

    return unpack(minimizer(result), range_nodes, range_wh, range_wl, D)
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

    return EmbeddedCubature(nodes, weights_high, weights_low)
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
