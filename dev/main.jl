using HAdaptiveIntegration: GrundmannMoeller, embedded_cubature
using Printf
using StaticArrays

function _gm_weight_next(weight::Vector{T}, dim::Int) where {T}
    s = length(weight)
    weight_next = Vector{T}(undef, s + 1)
    for i in s:-1:1
        v = dim + 1 + 2 * s - i
        weight_next[i + 1] = -weight[i] * (v - i)^2 / (4 * i * v)
    end
    v = dim + 1 + 2 * s
    weight_next[1] = -weight_next[2] * (T(v) / (v - 2))^(2 * s) / (v - 2)
    return weight_next
end

function multi_indexes(dim::Int, k_max::Int)
    if (dim ≤ 0) || (k_max < 0)
        return Vector{Vector{Tuple{}}}()
    end

    multi_indexes = [[tuple(k)] for k in 0:k_max]

    for d in 2:dim
        new = [Vector{NTuple{d,Int}}() for _ in 0:k_max]
        for (k, indexes) in zip(Iterators.countfrom(0), multi_indexes)
            for idx in indexes
                for n in 0:(k_max - k)
                    push!(new[k + 1 + n], (idx..., n))
                end
            end
        end
        multi_indexes = new
    end

    return multi_indexes
end

function embedded_cubature_new(T::DataType, gm::GrundmannMoeller{D}) where {D}
    oh, ol = gm.deg, gm.deg - 2
    sh, sl = (oh - 1) ÷ 2, (ol - 1) ÷ 2

    wl = [1 / T(factorial(D))]
    for _ in 1:sl
        wl = _gm_weight_next(wl, D)
    end
    wh = copy(wl)
    for _ in (sl + 1):sh
        wh = _gm_weight_next(wh, D)
    end

    nodes = Dict{NTuple{D,Rational{Int}},T}()
    for (i, indexes, w) in
        zip(Iterators.countfrom(0), multi_indexes(D + 1, sh), reverse(wh))
        dem = (D + 1 + 2 * i)
        for idx in indexes
            node = (2 .* idx[2:end] .+ 1) .// dem
            nodes[node] = w + get(nodes, node, 0//1)
        end
    end
    display(nodes)

    return nothing # EmbeddedCubature(nodes_high, weights_high, weights_low)
end

println("-"^64)
dim = 2
deg = 5

gm = GrundmannMoeller{dim}(deg)
ec = embedded_cubature(Rational{Int}, gm)
display(ec)

embedded_cubature_new(Rational{Int}, GrundmannMoeller{dim}(deg))
