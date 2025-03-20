using Base.Iterators: countfrom, partition
using BenchmarkTools
using ForwardDiff: jacobian
using HAdaptiveIntegration:
    EmbeddedCubature,
    SEGMENT_GK15,
    TETRAHEDRON_GM35,
    TRIANGLE_GM20,
    TRIANGLE_RL19,
    TabulatedEmbeddedCubature,
    embedded_cubature,
    integral_monomial_simplex
using LinearAlgebra
using StaticArrays

function pack(S::DataType, nodes, weights)
    U = Vector{S}()
    for (x, w) in zip(nodes, weights)
        append!(U, S.(x))
        push!(U, S(w))
    end
    return U
end

function unpack(U::Vector{T}, D::Int) where {T}
    nodes, weights = Vector{SVector{D,T}}(), Vector{T}()
    for t in partition(U, D + 1)
        push!(nodes, SVector{D}(t[1:D]))
        push!(weights, t[D + 1])
    end
    return nodes, weights
end

function increase_precision(
    S::DataType, ec::EmbeddedCubature{D,T}, order_high::Int, order_low::Int
) where {D,T}
    # n, w = increase_precision(S, ec.nodes, ec.weights_high, order_high)
    n, w = increase_precision(
        S, ec.nodes[1:length(ec.weights_low)], ec.weights_low, order_low
    )

    return n, w
end

function increase_precision(S::DataType, nodes, weights, order::Int)
    D = length(first(nodes))

    monomials, values = Vector{NTuple{D,Int}}(), Vector{S}()
    for pairs in integral_monomial_simplex(D, order)
        for (multi_index, value) in pairs
            push!(monomials, multi_index)
            push!(values, value)
        end
    end

    U = pack(S, nodes, weights)

    function F(u)
        v = zeros(eltype(u), length(monomials))
        for (i, mi) in enumerate(monomials)
            v[i] += sum(prod(t[1:D] .^ mi) * t[D + 1] for t in partition(u, D + 1))
        end
        v -= values
        return v
    end

    display(maximum(abs.(jacobian(F, U) \ F(U))))
    for _ in 1:10
        U -= jacobian(F, U) \ F(U)
        display(maximum(abs.(jacobian(F, U) \ F(U))))
    end

    println()
    display(U[1])
    display(abs(U[1] - 1 / 3))
    n, w = unpack(U, D)
    return n, w
end

# increase_precision(Float64, embedded_cubature(Float32, SEGMENT_GK15), 23, 13)

ec = embedded_cubature(Float64, TRIANGLE_RL19)
n, w = increase_precision(Float64, embedded_cubature(Float32, TRIANGLE_RL19), 7, 5)

println()
display(norm(norm.(n - ec.nodes[1:7], Inf), Inf))
display(norm(w - ec.weights_low, Inf))

# increase_precision(Float64, embedded_cubature(Float32, TETRAHEDRON_GM35), 7, 5)
