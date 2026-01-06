using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using Quadmath
using Test
using Base.Iterators: countfrom

function validate_orders(
    ec::EmbeddedCubature{D,T},
    ::Type{DOM},
    order_high::Int,
    order_low::Int;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T,DOM}
    val_ex = integral_monomials(DOM, order_high)
    for (k, α2v) in zip(countfrom(0), val_ex)
        for (α, v) in α2v
            vh, vl = eval_monomial(ec, α)

            vs = [(vh, "high order")]
            if k ≤ order_low
                push!(vs, (vl, "low order"))
            end

            for (v, label) in vs
                err_abs = abs(vh - v)
                err_rel = abs(err_abs / v)
                if err_abs > atol && err_rel > rtol
                    msg = "fail to integrate $label rule within tolerance at degree = $α.\n"
                    msg *= "Absolute error: $(err_abs), atol: $atol.\n"
                    msg *= "Relative error: $(err_rel), rtol: $rtol."
                    @error msg
                    return false
                end
            end
        end
    end

    return true
end

function integral_monomials(::Type{<:Segment}, deg_tot_max::Int)
    return [[(k,) => 1//(k + 1)] for k in 0:deg_tot_max]
end

function integral_monomials(::Type{<:Simplex{D}}, deg_tot_max::Int) where {D}
    @assert (D > 0) && (deg_tot_max ≥ 0) "must have `D > 0` and `k_max ≥ 0`."

    exponent2values = [[(i,) => 1//prod((i + 1):(i + D))] for i in 0:deg_tot_max]

    for d in 2:D
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, α2v) in zip(countfrom(0), exponent2values)
            for (α, v) in α2v
                push!(new[k + 1], (0, α...) => v)
                t = 1
                for n in 1:(deg_tot_max - k)
                    t *= n//(n + k + D)
                    push!(new[k + 1 + n], (n, α...) => t * v)
                end
            end
        end
        exponent2values = new
    end

    return exponent2values
end

function integral_monomials(::Type{<:Orthotope{D}}, deg_tot_max::Int) where {D}
    @assert (D > 0) && (deg_tot_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    exponent2values = [[(i,) => 1//(i + 1)] for i in 0:deg_tot_max]

    for d in 2:D
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, α2v) in zip(countfrom(0), exponent2values)
            for (α, v) in α2v
                for n in 0:(deg_tot_max - k)
                    push!(new[k + 1 + n], (n, α...) => v//(n + 1))
                end
            end
        end
        exponent2values = new
    end

    return exponent2values
end

function eval_monomial(ec::EmbeddedCubature{D,T}, α::NTuple{D,Int}) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)

    v = prod(ec.nodes[1] .^ α)
    Iₗ = ec.weights_low[1] * v
    Iₕ = ec.weights_high[1] * v
    for i in 2:L
        v = prod(ec.nodes[i] .^ α)
        Iₗ += ec.weights_low[i] * v
        Iₕ += ec.weights_high[i] * v
    end
    for i in (L + 1):H
        Iₕ += ec.weights_high[i] * prod(ec.nodes[i] .^ α)
    end

    return Iₕ, Iₗ
end

@testset "Rule construction" begin
    @testset "Tabulated embedded cubature" begin
        T = float(Int)

        tec = TabulatedEmbeddedCubature{Segment}(;
            description="Gauss (SEGMENT_G3)",
            reference="",
            precision=16,
            nodes=[
                [string(0.5)], [string((1 - √(3 / 5)) / 2)], [string((1 + √(3 / 5)) / 2)]
            ],
            weights_high=[string(4 / 9), string(5 / 18), string(5 / 18)],
            order_high=5,
            weights_low=[string(1.0)],
            order_low=1,
        )
        @test typeof(tec) <: TabulatedEmbeddedCubature{Segment}
        @test orders(tec) == (5, 1)

        ec = embedded_cubature(tec)
        @test typeof(ec) <: EmbeddedCubature{1,T}

        ec_ref = embedded_cubature(
            [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]],
            [4 / 9, 5 / 18, 5 / 18],
            [1.0],
        )
        @test typeof(ec_ref) <: EmbeddedCubature{1,T}

        @test ec.nodes ≈ ec_ref.nodes
        @test ec.weights_high ≈ ec_ref.weights_high
        @test ec.weights_low ≈ ec_ref.weights_low
    end

    @testset "Radon-Laurie" begin
        rl = RadonLaurie()
        @test typeof(rl) <: AbstractRule{Simplex{2}}

        @test orders(rl) == (8, 5)
        @test validate_orders(embedded_cubature(rl), Triangle, 8, 5)
    end

    @testset "Grundmann-Möller" begin
        @test_throws AssertionError typeof(GrundmannMoeller{1}(7, 6))
        @test_throws AssertionError typeof(GrundmannMoeller{1}(5, 7))

        gm = GrundmannMoeller{2}(5, 3)
        @test typeof(gm) <: AbstractRule{Simplex{2}}
        @test orders(gm) == (5, 3)
        @test validate_orders(embedded_cubature(gm), Triangle, 5, 3)
    end

    @testset "Genz-Malik" begin
        gm = GenzMalik{2}()
        @test typeof(gm) <: AbstractRule{Orthotope{2}}
        @test orders(gm) == (7, 5)
        @test validate_orders(embedded_cubature(gm), Rectangle, 7, 5)
    end
end

@testset "Tabulated rules" begin
    T = Float128

    @testset "Segment" begin
        for tec in (SEGMENT_GK7, SEGMENT_GK15, SEGMENT_GK31)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, Segment, orders(tec)...)
        end

        ec = embedded_cubature(GrundmannMoeller{1}(7, 5), T)
        @test validate_orders(ec, Simplex{1}, 7, 5)

        ec = embedded_cubature(GenzMalik{1}(), T)
        @test validate_orders(ec, Orthotope{1}, 7, 5)
    end

    @testset "Square" begin
        for tec in (SQUARE_CH21, SQUARE_CH25)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, Rectangle, orders(tec)...)
        end
    end

    @testset "Cube" begin
        for tec in (CUBE_BE65, CUBE_BE115)
            ec = embedded_cubature(tec, T)
            @test validate_orders(ec, Cuboid, orders(tec)...)
        end
    end
end
