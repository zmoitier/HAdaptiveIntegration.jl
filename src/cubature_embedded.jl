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

function embedded_cubature(
    T::DataType, tec::TabulatedEmbeddedCubature{DOM}
) where {DOM<:AbstractDomain}
    if 10 * eps(T) < 10.0^(-tec.nb_significant_digits)
        @warn "The embedded cubature `$(tec.description)` has fewer significant digits than type $T, which may lead to numerical inaccuracies in computations."
    end
    D = dimension(DOM)
    return EmbeddedCubature(
        [SVector{D,T}(parse.(T, x)) for x in tec.nodes],
        parse.(T, tec.weights_high),
        parse.(T, tec.weights_low),
    )
end
embedded_cubature(tec::TabulatedEmbeddedCubature) = embedded_cubature(float(Int), tec)

# Based on Eq. (4.1) in the article:
#   A. Grundmann and H. M. Möller, Invariant integration formulas for the n-simplex by
#   combinatorial methods, SIAM Journal on Numerical Analysis, volume 15, 1978,
#   https://doi.org/10.1137/0715019.
function embedded_cubature(T::DataType, gm::GrundmannMoeller{D}) where {D}
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
embedded_cubature(gm::GrundmannMoeller) = embedded_cubature(float(Int), gm)

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

# Based on the article:
#   A. C. Genz and A. A. Malik, Remarks on algorithm 006: An adaptive algorithm for
#   numerical integration over an N-dimensional rectangular region, Journal of Computational
#   and Applied Mathematics, Volume 6, Issue 4, 1980,
#   https://doi.org/10.1016/0771-050X(80)90039-X.
function embedded_cubature(T::DataType, ::GenzMalik{D}) where {D}
    # λ₁                 # (0, ..., 0)
    λ₂ = sqrt(T(9) / 70) # (λ₂, 0, ..., 0)
    λ₃ = sqrt(T(9) / 10) # (λ₃, 0, ..., 0)
    λ₄ = λ₃              # (λ₄, λ₄, 0, ..., 0)
    λ₅ = sqrt(T(9) / 19) # (λ₅, ..., λ₅)

    wh = SVector{5,T}(
        (12_824 - 9_120 * D + 400 * D^2)//19_683,
        980//6_561,
        (1_820 - 400 * D)//19_683,
        200//19_683,
        6_859//(19_683 * 2^D),
    )
    # w' weights
    wl = SVector{4,T}(
        (729 - 950 * D + 50 * D^2)//729, 245//486, (265 - 100 * D)//1_458, 25//729
    )

    nodes = [zero(SVector{D,T})]
    weights_high = [wh[1]]
    weights_low = [wl[1]]

    node = zero(MVector{D,T})
    for i in 1:D
        for (λ, wₕ, wₗ) in zip((λ₂, λ₃), wh[2:3], wl[2:3])
            for s in (1, -1)
                node[i] = s * λ
                push!(nodes, SVector{D,T}(node))
                push!(weights_high, wₕ)
                push!(weights_low, wₗ)
            end
        end
        node[i] = 0
    end

    for i in 1:(D - 1)
        for j in (i + 1):D
            for (s₁, s₂) in Iterators.product((1, -1), (1, -1))
                node[i] = s₁ * λ₄
                node[j] = s₂ * λ₄
                push!(nodes, SVector{D,T}(node))
                push!(weights_high, wh[4])
                push!(weights_low, wl[4])
            end
            node[j] = 0
        end
        node[i] = 0
    end

    node .= λ₅
    for signs in Iterators.product([(1, -1) for _ in 1:D]...)
        push!(nodes, SVector{D,T}(signs .* node))
        push!(weights_high, wh[5])
    end

    Φ = x -> (x .+ 1) ./ 2
    return EmbeddedCubature(Φ.(nodes), weights_high, weights_low)
end
embedded_cubature(gm::GenzMalik{D}) where {D} = embedded_cubature(float(Int), gm)

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
    fct, domain::AbstractDomain{D,T}, norm=x -> LinearAlgebra.norm(x, Inf)
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
