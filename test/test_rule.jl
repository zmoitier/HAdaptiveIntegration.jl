using Test

import HAdaptiveIntegration as hai

@testset "Abstract rule" begin
    @testset "Tabulated embedded cubature" begin
        tec = hai.TabulatedEmbeddedCubature{hai.Segment}(;
            description="Gauss (SEGMENT_G3)",
            reference="",
            nb_significant_digits=16,
            nodes=[["5e-1"], ["1.127016653792583e-1"], ["8.872983346207417e-1"]],
            weights_high=[
                "4.444444444444444e-1", "2.777777777777778e-1", "2.777777777777778e-1"
            ],
            order_high=5,
            weights_low=["1"],
            order_low=1,
        )
        @test typeof(tec) <: hai.TabulatedEmbeddedCubature{hai.Segment}
        @test hai.orders(tec) == (5, 1)

        ec = hai.embedded_cubature(tec)
        @test typeof(ec) <: hai.EmbeddedCubature{1,Float64}

        ec_ref = hai.embedded_cubature(
            [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]],
            [4 / 9, 5 / 18, 5 / 18],
            [1.0],
        )
        @test typeof(ec_ref) <: hai.EmbeddedCubature{1,Float64}

        @test ec.nodes ≈ ec_ref.nodes
        @test ec.weights_high ≈ ec_ref.weights_high
        @test ec.weights_low ≈ ec_ref.weights_low
    end

    @testset "Grundmann-Möller" begin
        @test typeof(hai.GrundmannMoeller{3}(7, 3)) <: hai.AbstractRule{hai.Simplex{3}}
        @test_throws AssertionError typeof(hai.GrundmannMoeller{3}(7, 6))
        @test_throws AssertionError typeof(hai.GrundmannMoeller{3}(5, 7))

        gm = hai.GrundmannMoeller{4}(5, 3)
        @test hai.orders(gm) == (5, 3)
        @test hai.validate_orders(
            hai.embedded_cubature(Float64, gm),
            hai.Simplex,
            gm.order_high,
            gm.order_low;
            rtol=20 * eps(Float64),
        )
    end

    @testset "Genz-Malik" begin
        @test typeof(hai.GenzMalik{2}()) <: hai.AbstractRule{hai.Orthotope{2}}
        @test typeof(hai.GenzMalik{3}()) <: hai.AbstractRule{hai.Orthotope{3}}

        gm = hai.GenzMalik{4}()
        @test hai.validate_orders(
            hai.embedded_cubature(Float64, gm),
            hai.Orthotope,
            hai.orders(gm)...;
            rtol=10 * eps(Float64),
        )
    end
end

@testset "Tabulated rules" begin
    T = Float64

    @testset "Segment" begin
        for tec in (
            # hai.SEGMENT_GK7,
            hai.SEGMENT_GK15,
            hai.SEGMENT_GK31,
        )
            ec = hai.embedded_cubature(T, tec)
            @test hai.validate_orders(ec, hai.Orthotope, hai.orders(tec)...)
        end
    end

    @testset "Square" begin
        for tec in (hai.SQUARE_CH25, hai.SQUARE_CH21, hai.SQUARE_GM17)
            ec = hai.embedded_cubature(T, tec)
            @test hai.validate_orders(ec, hai.Orthotope, hai.orders(tec)...)
        end
    end

    @testset "Triangle" begin
        for tec in (hai.TRIANGLE_RL19, hai.TRIANGLE_GM19)
            ec = hai.embedded_cubature(T, tec)
            @test hai.validate_orders(ec, hai.Simplex, hai.orders(tec)...)
        end
    end

    @testset "Cube" begin
        for tec in (hai.CUBE_BE65, hai.CUBE_GM33)
            ec = hai.embedded_cubature(T, tec)
            @test hai.validate_orders(ec, hai.Orthotope, hai.orders(tec)...)
        end
    end

    @testset "Tetrahedron" begin
        tec = hai.TETRAHEDRON_GM35
        ec = hai.embedded_cubature(T, tec)
        @test hai.validate_orders(ec, hai.Simplex, hai.orders(tec)..., rtol=20 * eps(T))
    end
end
