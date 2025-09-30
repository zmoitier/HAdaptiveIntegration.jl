"""
    struct RadonLaurie <: AbstractRule{Simplex{2}}

Embedded cubature rule for a `2`-simplex (*a.k.a* a triangle) of high order `8` and low
order `5`.
"""
struct RadonLaurie <: AbstractRule{Simplex{2}} end

function orders(::RadonLaurie)
    return 8, 5
end

# Based on the article:
#   D. P. Laurie. 1982. Algorithm 584: CUBTRI: Automatic Cubature over a Triangle. ACM
#   Trans. Math. Softw. 8, 1982, https://doi.org/10.1145/355993.356001.
function embedded_cubature(::RadonLaurie, (::Type{T})=float(Int)) where {T}
    ϕ, σ = sqrt(T(15)), sqrt(T(7))

    nodes = [SVector{2,T}(fill(1//3, 2))]
    weights_high = T[(7_137 // 62_720 - σ * 45 // 1568) / 2]
    weights_low = T[9 // 40 / 2]

    for s in (1, -1)
        ζ₁ = 3//7 + s * ϕ * 2//21
        ζ₂ = 2//7 - s * ϕ * 1//21
        wₕ =
            (
                -9_301_697//469_5040 - s * ϕ * 13_517_313//23_475_200 +
                σ * 764_885//939_008 +
                s * ϕ * σ * 198_763//939_008
            ) / 6
        wₗ = (31//80 - s * ϕ * 1//400) / 6

        append!(nodes, SVector{2,T}.([(ζ₁, ζ₂), (ζ₂, ζ₁), (ζ₂, ζ₂)]))
        append!(weights_high, [wₕ, wₕ, wₕ])
        append!(weights_low, [wₗ, wₗ, wₗ])
    end

    for s in (1, -1)
        ζ₁ = 4//9 + s * ϕ * 1//9 + σ * 1//9 - s * ϕ * σ * 1//45
        ζ₂ = 5//18 - s * ϕ * 1//18 - σ * 1//18 + s * ϕ * σ * 1//90
        wₕ =
            (
                102_791_225//5_915_7504 + s * ϕ * 23_876_225//59_157_504 -
                σ * 34_500_875//59_157_504 - s * ϕ * σ * 9_914_825//59_157_504
            ) / 6

        append!(nodes, SVector{2,T}.([(ζ₁, ζ₂), (ζ₂, ζ₁), (ζ₂, ζ₂)]))
        append!(weights_high, [wₕ, wₕ, wₕ])
    end

    ζ₁ = 5//18 - ϕ * 1//18 - σ * 1//18 + ϕ * σ * 1//90
    ζ₂ = 5//18 + ϕ * 1//18 - σ * 1//18 - ϕ * σ * 1//90
    ζ₃ = 4//9 + σ / 9
    wₕ = (11_075//8_064 - σ * 125//288) / 12
    append!(
        nodes, SVector{2,T}.([(ζ₁, ζ₂), (ζ₁, ζ₃), (ζ₂, ζ₁), (ζ₂, ζ₃), (ζ₃, ζ₁), (ζ₃, ζ₂)])
    )
    append!(weights_high, [wₕ, wₕ, wₕ, wₕ, wₕ, wₕ])

    return EmbeddedCubature(nodes, weights_high, weights_low)
end
