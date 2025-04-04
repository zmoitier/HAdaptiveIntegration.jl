using Base.Iterators: countfrom, partition
using BenchmarkTools
using ForwardDiff: jacobian
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
using Optim
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

# [ nodes[1]..., nodes[2]..., ..., nodes[end]..., weights_high..., weights_low... ]
function pack(S::DataType, ec::EmbeddedCubature{D,T}) where {D,T}
    U = Vector{S}()
    for node in ec.nodes
        append!(U, node)
    end
    append!(U, ec.weights_high)
    append!(U, ec.weights_low)
    return U
end

function unpack(U::Vector{T}, D::Int, H::Int) where {T}
    nodes = [SVector{D}(t) for t in partition(U[1:(H * D)], D)]
    weights_high = U[(H * D + 1):(H * D + H)]
    weights_low = U[(H * D + H + 1):end]
    return EmbeddedCubature(nodes, weights_high, weights_low)
end

function increase_precision(
    S::DataType, ec::EmbeddedCubature{D,T}, order_high::Int, order_low::Int
) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)

    monomials, values = Vector{NTuple{D,Int}}(), Vector{S}()
    for pairs in integral_chebyshev_orthotope(D, order_high)
        for (multi_index, value) in pairs
            push!(monomials, multi_index)
            push!(values, value)
        end
    end

    function F(u)
        v = zero(eltype(u))

        nodes = @view u[1:(H * D)]
        weights_hi = @view u[(H * D + 1):(H * D + H)]
        weights_lo = @view u[(H * D + H + 1):end]

        acos_nodes = acos.(2 .* nodes .- 1)

        for (mi, val) in zip(monomials[1:order_low], values[1:order_low])
            v_lo, v_hi = zero(v), zero(v)
            for (i, acn) in enumerate(partition(acos_nodes, D))
                tmp = prod(cos.(mi .* acn))
                v_hi += weights_hi[i] * tmp
                if i ≤ L
                    v_lo += weights_lo[i] * tmp
                end
            end
            v += (v_hi - val)^2 + (v_lo - val)^2
        end
        for (mi, val) in zip(monomials[(order_low + 1):end], values[(order_low + 1):end])
            v_hi = zero(v)
            for (acn, w_hi) in zip(partition(acos_nodes, D), weights_hi)
                v_hi += w_hi * prod(cos.(mi .* acn))
            end
            v += (v_hi - val)^2
        end

        return v
    end

    # function F(u)
    #     v = zero(eltype(u))

    #     nodes = @view u[1:(H * D)]
    #     weights_hi = @view u[(H * D + 1):(H * D + H)]
    #     weights_lo = @view u[(H * D + H + 1):end]

    #     for (mi, val) in zip(monomials[1:order_low], values[1:order_low])
    #         v_lo, v_hi = zero(v), zero(v)
    #         for (i, n) in enumerate(partition(nodes, D))
    #             tmp = prod(n .^ mi)
    #             v_hi += weights_hi[i] * tmp
    #             if i ≤ L
    #                 v_lo += weights_lo[i] * tmp
    #             end
    #         end
    #         v += (v_hi - val)^2 + (v_lo - val)^2
    #     end
    #     for (mi, val) in zip(monomials[(order_low + 1):end], values[(order_low + 1):end])
    #         v_hi = zero(v)
    #         for (n, w_hi) in zip(partition(nodes, D), weights_hi)
    #             v_hi += w_hi * prod(n .^ mi)
    #         end
    #         v += (v_hi - val)^2
    #     end

    #     return v
    # end

    U = pack(S, ec)
    U += randn(S, length(U)) * 1e-8
    @show F(U)

    @show result = optimize(
        F, U, LBFGS(), Optim.Options(; x_tol=1e-40, g_tol=1e-50); autodiff=:forward
    )
    V = Optim.minimizer(result)

    return unpack(V, D, length(ec.weights_high))
end

println("-"^32)

# tec = SEGMENT_GK15
# tec = SEGMENT_GK31
# tec = SQUARE_CHG21
tec = SQUARE_CHG25
# tec = TRIANGLE_RL19
# tec = CUBE_BE65
# tec = TETRAHEDRON_GM35

H, L = length(tec.weights_high), length(tec.weights_low)
ec = embedded_cubature(Float64, tec)
ec_ex = increase_precision(
    BigFloat, embedded_cubature(Float64, tec), tec.order_high, tec.order_low + 2
)

display(norm(norm.(ec_ex.nodes - ec.nodes, Inf), Inf))
println("Low order")
display(norm(ec_ex.weights_low - ec.weights_low, Inf))
println()
println("High order")
display(norm(ec_ex.weights_high - ec.weights_high, Inf))

display(ec_ex.weights_low)
