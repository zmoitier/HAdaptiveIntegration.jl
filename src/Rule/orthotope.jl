"""
   struct GenzMalik{D} <: AbstractRule{Orthotope{D}}

Embedded cubature rule for a `D`-orthotope of high order `7` and low order `5`.

## Type Parameters:
- `D`: The dimension of the orthotope.
"""
struct GenzMalik{D} <: AbstractRule{Orthotope{D}} end

function orders(::GenzMalik)
    return 7, 5
end

# Based on the article:
#   A. C. Genz and A. A. Malik, Remarks on algorithm 006: An adaptive algorithm for
#   numerical integration over an N-dimensional rectangular region, Journal of Computational
#   and Applied Mathematics, Volume 6, Issue 4, 1980,
#   https://doi.org/10.1016/0771-050X(80)90039-X.
function embedded_cubature(T::DataType, ::GenzMalik{D}) where {D}
    # map to the reference domain
    Φ = x -> (x .+ 1) ./ 2

    # w₁, ..., w₅ for the reference domain (divided by 2^D)
    wh = SVector{5,T}(
        (12_824 - 9_120 * D + 400 * D^2)//19_683,
        980//6_561,
        (1_820 - 400 * D)//19_683,
        200//19_683,
        6_859//(19_683 * 2^D),
    )
    # w'₁, ..., w'₄ for the reference domain (divided by 2^D)
    wl = SVector{4,T}(
        (729 - 950 * D + 50 * D^2)//729, 245//486, (265 - 100 * D)//1_458, 25//729
    )

    # start with the node (0, ..., 0)
    node = zeros(T, D)
    nodes = [Φ(SVector{D,T}(node))]
    weights_high = [wh[1]]
    weights_low = [wl[1]]

    # generate orbit of the point (λ₂, 0, ..., 0) and (λ₃, 0, ..., 0)
    λ₂ = sqrt(T(9) / 70)
    λ₃ = sqrt(T(9) / 10)
    for (λ, wₕ, wₗ) in zip((λ₂, λ₃), wh[2:3], wl[2:3])
        for i in 1:D
            for s in (1, -1)
                node[i] = s * λ
                push!(nodes, Φ(SVector{D,T}(node)))
                push!(weights_high, wₕ)
                push!(weights_low, wₗ)
            end
            node[i] = 0
        end
    end

    # generate orbit of the point (λ₄, λ₄, 0, ..., 0)
    λ₄ = λ₃
    for i in 1:(D - 1)
        for j in (i + 1):D
            for (s₁, s₂) in Iterators.product((1, -1), (1, -1))
                node[i] = s₁ * λ₄
                node[j] = s₂ * λ₄
                push!(nodes, Φ(SVector{D,T}(node)))
                push!(weights_high, wh[4])
                push!(weights_low, wl[4])
            end
            node[j] = 0
        end
        node[i] = 0
    end

    # generate orbit of the point (λ₅, λ₅, ..., λ₅)
    λ₅ = sqrt(T(9) / 19)
    node .= λ₅
    for signs in Iterators.product([(1, -1) for _ in 1:D]...)
        push!(nodes, Φ(SVector{D,T}(signs .* node)))
        push!(weights_high, wh[5])
    end

    return EmbeddedCubature(nodes, weights_high, weights_low)
end
embedded_cubature(gm::GenzMalik{D}) where {D} = embedded_cubature(float(Int), gm)
