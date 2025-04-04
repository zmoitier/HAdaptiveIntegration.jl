using Base.Iterators: countfrom, partition
using BenchmarkTools
using ForwardDiff: jacobian
using GenericLinearAlgebra
using HAdaptiveIntegration:
    CUBE_BE65,
    EmbeddedCubature,
    SEGMENT_GK15,
    SEGMENT_GK31,
    SQUARE_CHG21,
    SQUARE_CHG25,
    TETRAHEDRON_GM35,
    TRIANGLE_GM20,
    TRIANGLE_RL19,
    TabulatedEmbeddedCubature,
    embedded_cubature,
    integral_monomial_orthotope,
    integral_monomial_simplex
using LinearAlgebra
using SparseArrays
using StaticArrays

function integral_chebyshev_orthotope(d::Int, tdm::Int)
    if d ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(k,) => iseven(k) ? 1//(1 - k * k) : 0//1] for k in 0:tdm]

    for n in 2:d
        tmp = [Vector{Pair{NTuple{n,Int},Rational{Int}}}() for _ in 0:tdm]
        for (td, idx_val) in zip(Iterators.countfrom(0), indexes)
            for (idx, val) in idx_val
                for k in 0:(tdm - td)
                    push!(
                        tmp[td + 1 + k],
                        (k, idx...) => val * (iseven(k) ? 1//(1 - k * k) : 0//1),
                    )
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

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
    nlo, wlo = increase_precision(
        S, ec.nodes[1:length(ec.weights_low)], ec.weights_low, order_low
    )
    nhi, whi = increase_precision(S, ec.nodes, ec.weights_high, order_high)

    return nlo, wlo, nhi, whi
end

function increase_precision(S::DataType, nodes, weights, order::Int)
    D = length(first(nodes))

    monomials, values = Vector{NTuple{D,Int}}(), Vector{S}()
    for pairs in integral_chebyshev_orthotope(D, order)
        for (multi_index, value) in pairs
            push!(monomials, multi_index)
            push!(values, value)
        end
    end

    U = pack(S, nodes, weights)

    function F(u)
        v = zeros(eltype(u), length(monomials))
        for (i, mi) in enumerate(monomials)
            v[i] += sum(
                prod(cos.(mi .* acos.(2 .* t[1:D] .- 1))) * t[D + 1] for
                t in partition(u, D + 1)
            )
        end
        v -= values
        return v
    end

    f, J = F(U), jacobian(F, U)

    qr_decomp = qr(J)
    display(qr_decomp.Q)
    display(sparse(qr_decomp.R))

    @show norm(f, Inf)
    @show size(J)
    @show cond(J)

    println("start newton")
    display(maximum(abs.(jacobian(F, U) \ F(U))))
    for _ in 1:3
        U -= jacobian(F, U) \ F(U)
        display(maximum(abs.(jacobian(F, U) \ F(U))))
    end
    println("end newton")

    println()

    n, w = unpack(U, D)
    return n, w
end

println("-"^32)

tec = SEGMENT_GK15
# tec = SEGMENT_GK31
# tec = SQUARE_CHG21
# tec = SQUARE_CHG25
# tec = TRIANGLE_RL19
# tec = CUBE_BE65
# tec = TETRAHEDRON_GM35

H, L = length(tec.weights_high), length(tec.weights_low)
ec = embedded_cubature(Float64, tec)
nlo, wlo, nhi, whi = increase_precision(
    Float64, embedded_cubature(Float32, tec), tec.order_high, tec.order_low
)

println("Low order")
display(norm(norm.(nlo - ec.nodes[1:L], Inf), Inf))
display(norm(wlo - ec.weights_low, Inf))
println()
println("High order")
display(norm(norm.(nhi - ec.nodes, Inf), Inf))
display(norm(whi - ec.weights_high, Inf))
