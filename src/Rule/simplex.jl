"""
    struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}

Embedded cubature rule for a `D`-simplex.

## Type Parameters:
- `D`: The dimension of the simplex.

## Fields:
- `order_high::Int`: the high order of the rule.
- `order_low::Int`: the low order of the rule.

## Invariants (check at construction):
- `order_high` and `order_low` must be odd.
- must have `order_high > order_low ≥ 1`.
"""
struct GrundmannMoeller{D} <: AbstractRule{Simplex{D}}
    order_high::Int
    order_low::Int

    function GrundmannMoeller{D}(order_high::Int, order_low::Int) where {D}
        @assert isodd(order_high) && isodd(order_low)
        @assert order_high > order_low ≥ 1
        return new{D}(order_high, order_low)
    end
end

function orders(gm::GrundmannMoeller)
    return gm.order_high, gm.order_low
end

# Based on Eq. (4.1) in the article:
#   A. Grundmann and H. M. Möller, Invariant integration formulas for the n-simplex by
#   combinatorial methods, SIAM Journal on Numerical Analysis, volume 15, 1978,
#   https://doi.org/10.1137/0715019.
function embedded_cubature(gm::GrundmannMoeller{D}, (::Type{T})=float(Int)) where {D,T}
    sh, sl = (gm.order_high - 1) ÷ 2, (gm.order_low - 1) ÷ 2

    # Grundmann-Möller weights computed iteratively instead of using the formula for
    # numerical stability at higher degree.
    gm_wl = [1 / reduce(*, 1:D; init=T(1))]
    for _ in 1:sl
        gm_wl = _gm_weight_next(gm_wl, D)
    end
    gm_wh = copy(gm_wl)
    for _ in (sl + 1):sh
        gm_wh = _gm_weight_next(gm_wh, D)
    end

    nodes = Vector{SVector{D,T}}()
    weights_high = Vector{T}()
    weights_low = Vector{T}()

    count = 0
    node_to_idx = Dict{NTuple{D,Rational{Int}},Int}()

    mlt_idx = _multi_indexes(D + 1, sh)
    for i in 0:sh
        dem = D + 1 + 2 * i
        for idx in mlt_idx[i + 1]
            node = (2 .* idx[2:end] .+ 1) .// dem

            if node in keys(node_to_idx)
                j = node_to_idx[node]
                weights_high[j] += gm_wh[sh - i + 1]
                if i ≤ sl
                    weights_low[j] += gm_wl[sl - i + 1]
                end
            else
                count += 1
                node_to_idx[node] = count

                push!(nodes, SVector{D,T}(node))
                push!(weights_high, gm_wh[sh - i + 1])
                if i ≤ sl
                    push!(weights_low, gm_wl[sl - i + 1])
                end
            end
        end
    end

    return EmbeddedCubature(nodes, weights_high, weights_low)
end

function _gm_weight_next(weight::Vector{T}, dim::Int) where {T}
    s = length(weight)

    weight_next = Vector{T}(undef, s + 1)
    for i in s:-1:1
        v = dim + 1 + 2 * s - i
        weight_next[i + 1] = -weight[i] * (v - i)^2 / (4 * i * v)
    end
    v = dim + 2 * s - 1
    weight_next[1] = -weight_next[2] * (T(v + 2) / v)^(2 * s) / v

    return weight_next
end

function _multi_indexes(dim::Int, k_max::Int)
    @assert (dim > 0) && (k_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    mlt_idx = [[tuple(k)] for k in 0:k_max]

    for d in 2:dim
        new = [Vector{NTuple{d,Int}}() for _ in 0:k_max]
        for (k, mis) in zip(Iterators.countfrom(0), mlt_idx)
            for mi in mis
                for n in 0:(k_max - k)
                    push!(new[k + 1 + n], (mi..., n))
                end
            end
        end
        mlt_idx = new
    end

    return mlt_idx
end
