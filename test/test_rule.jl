using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Polynomial
using HAdaptiveIntegration.Rule
using Quadmath
using Test

function validate_orders(
    ec::EmbeddedCubature{D,T},
    ::Type{DOM},
    order_high::Int,
    order_low::Int;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T,DOM}
    val_ex = integral_monomials_exact(DOM, order_high)
    val_num_hi, val_num_lo = integral_monomials_rule(ec, order_high, order_low)

    for k in 0:order_high
        for (α2v_num, α2v_ex) in zip(val_num_hi[k + 1], val_ex[k + 1])
            for ((α, v_num), (_, v_ex)) in zip(α2v_num, α2v_ex)
                err_abs = abs(v_num - v_ex)
                err_rel = abs(err_abs / v_ex)
                if err_abs > atol && err_rel > rtol
                    msg = "fail to integrate within tolerance at degree = $α.\n"
                    msg *= "Absolute error: $(err_abs), atol: $atol.\n"
                    msg *= "Relative error: $(err_rel), rtol: $rtol."
                    @error msg
                    return false
                end
            end
        end
    end

    for k in 0:order_low
        for (α2v_num, α2v_ex) in zip(val_num_lo[k + 1], val_ex[k + 1])
            for ((α, v_num), (_, v_ex)) in zip(α2v_num, α2v_ex)
                err_abs = abs(v_num - v_ex)
                err_rel = abs(err_abs / v_ex)
                if err_abs > atol && err_rel > rtol
                    msg = "fail to integrate within tolerance at degree = $α.\n"
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
