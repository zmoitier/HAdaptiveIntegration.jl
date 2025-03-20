using BenchmarkTools
using HAdaptiveIntegration:
    EmbeddedCubature,
    SEGMENT_GK15,
    TETRAHEDRON_GM35,
    TRIANGLE_GM20,
    TabulatedEmbeddedCubature,
    embedded_cubature
using LinearAlgebra
using StaticArrays

function integral_monomial_simplex(d::Int, tdm::Int)
    if d ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(i,) => 1//prod((i + 1):(i + d))] for i in 0:tdm]

    for n in 2:d
        tmp = [Vector{Pair{NTuple{n,Int},Rational{Int}}}() for _ in 0:tdm]
        for (td, idx_val) in zip(Iterators.countfrom(0), indexes)
            for (idx, val) in idx_val
                push!(tmp[td + 1], (0, idx...) => val)
                v = 1
                for k in 1:(tdm - td)
                    v *= k//(k + td + d)
                    push!(tmp[td + 1 + k], (k, idx...) => val * v)
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

function pack(
    S::DataType, nodes::SVector{N,SVector{D,T}}, weights::SVector{N,T}
) where {N,D,T}
    # U = [x₁..., ..., xₙ..., w...] where n = D
    U = Vector{S}(undef, (D + 1) * N)
    for (i, (x, w)) in enumerate(zip(nodes, weights))
        for n in 1:D
            U[(n - 1) * N + i] = S(x[n])
        end
        U[D * N + i] = S(w)
    end
    return U
end

function unpack(U::Vector{T}, D::Int) where {T}
    N = length(U) ÷ (D + 1)
    nodes, weights = Vector{SVector{D,T}}(), Vector{T}()
    for i in 1:N
        push!(nodes, SVector{D}(U[(n - 1) * N + i] for n in 1:D))
        push!(weights, U[D * N + i])
    end
    return nodes, weights
end

function increase_precision(
    S::DataType, ec::EmbeddedCubature{H,L,D,T}, order_high::Int, order_low::Int
) where {H,L,D,T}
    increase_precision(S, ec.nodes, ec.weights_high, order_high)
    # increase_precision(S, ec.nodes[1:L], ec.weights_low, order_low)

    return nothing
end

function increase_precision(
    S::DataType, nodes::SVector{N,SVector{D,T}}, weights::SVector{N,T}, order::Int
) where {N,D,T}
    monomials, values = Vector{NTuple{D,Int}}(), Vector{S}()
    for pairs in integral_monomial_simplex(D, order)
        for (multi_index, value) in pairs
            push!(monomials, multi_index)
            push!(values, value)
        end
    end

    U = pack(S, nodes, weights)

    # n, w = unpack(U, D)
    return nothing
end

# increase_precision(Float64, embedded_cubature(Float32, SEGMENT_GK15), 23, 13)
increase_precision(Float64, embedded_cubature(Float32, TRIANGLE_GM20), 7, 5)
# increase_precision(Float64, embedded_cubature(Float32, TETRAHEDRON_GM35), 7, 5)
