using HAdaptiveIntegration:
    AbstractRule,
    CUBE_BE115,
    CUBE_BE65,
    CUBE_GM33,
    EmbeddedCubature,
    GenzMalik,
    GrundmannMoeller,
    Orthotope,
    RadonLaurie,
    SEGMENT_GK15,
    SEGMENT_GK31,
    SEGMENT_GK7,
    SQUARE_CH21,
    SQUARE_CH25,
    SQUARE_GM17,
    Segment,
    Simplex,
    TETRAHEDRON_GM35,
    TRIANGLE_GM19,
    TRIANGLE_RL19,
    TabulatedEmbeddedCubature,
    embedded_cubature,
    orders,
    validate_orders
using Quadmath
using Test

@testset "Abstract rule" begin
    @testset "Tabulated embedded cubature" begin
        tec = TabulatedEmbeddedCubature{Segment}(;
            description="Gauss (SEGMENT_G3)",
            reference="",
            precision=16,
            nodes=[["5e-1"], ["1.127016653792583e-1"], ["8.872983346207417e-1"]],
            weights_high=[
                "4.444444444444444e-1", "2.777777777777778e-1", "2.777777777777778e-1"
            ],
            order_high=5,
            weights_low=["1"],
            order_low=1,
        )
        @test typeof(tec) <: TabulatedEmbeddedCubature{Segment}
        @test orders(tec) == (5, 1)

        ec = embedded_cubature(tec)
        @test typeof(ec) <: EmbeddedCubature{1,Float64}

        ec_ref = embedded_cubature(
            [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]],
            [4 / 9, 5 / 18, 5 / 18],
            [1.0],
        )
        @test typeof(ec_ref) <: EmbeddedCubature{1,Float64}

        @test ec.nodes ≈ ec_ref.nodes
        @test ec.weights_high ≈ ec_ref.weights_high
        @test ec.weights_low ≈ ec_ref.weights_low
    end

    @testset "Radon-Laurie" begin
        rl = RadonLaurie()
        @test typeof(rl) <: AbstractRule{Simplex{2}}

        @test orders(rl) == (8, 5)
        @test validate_orders(
            embedded_cubature(rl), Simplex, orders(rl)...; rtol=10 * eps(float(Int))
        )
    end

    @testset "Grundmann-Möller" begin
        @test typeof(GrundmannMoeller{3}(7, 3)) <: AbstractRule{Simplex{3}}
        @test_throws AssertionError typeof(GrundmannMoeller{3}(7, 6))
        @test_throws AssertionError typeof(GrundmannMoeller{3}(5, 7))

        gm = GrundmannMoeller{4}(5, 3)
        @test orders(gm) == (5, 3)
        @test validate_orders(
            embedded_cubature(gm), Simplex, orders(gm)...; rtol=20 * eps(float(Int))
        )
    end

    @testset "Genz-Malik" begin
        @test typeof(GenzMalik{2}()) <: AbstractRule{Orthotope{2}}
        @test typeof(GenzMalik{3}()) <: AbstractRule{Orthotope{3}}

        gm = GenzMalik{4}()
        @test validate_orders(
            embedded_cubature(gm), Orthotope, orders(gm)...; rtol=10 * eps(float(Int))
        )
    end
end

@testset "Tabulated rules" begin
    T = Float128
    rtol = 10 * eps(T)

    @testset "Segment" begin
        for tec in (SEGMENT_GK7, SEGMENT_GK15, SEGMENT_GK31)
            ec = embedded_cubature(T, tec)
            @test validate_orders(ec, Orthotope, orders(tec)...; rtol=rtol)
        end
    end

    @testset "Triangle" begin
        for tec in (TRIANGLE_GM19, TRIANGLE_RL19)
            ec = embedded_cubature(T, tec)
            @test validate_orders(ec, Simplex, orders(tec)...; rtol=rtol)
        end
    end

    @testset "Square" begin
        for tec in (SQUARE_GM17, SQUARE_CH21, SQUARE_CH25)
            ec = embedded_cubature(T, tec)
            @test validate_orders(ec, Orthotope, orders(tec)...; rtol=rtol)
        end
    end

    @testset "Tetrahedron" begin
        tec = TETRAHEDRON_GM35
        ec = embedded_cubature(T, tec)
        @test validate_orders(ec, Simplex, orders(tec)...; rtol=rtol)
    end

    @testset "Cube" begin
        for tec in (CUBE_GM33, CUBE_BE65, CUBE_BE115)
            ec = embedded_cubature(T, tec)
            @test validate_orders(ec, Orthotope, orders(tec)...; rtol=rtol)
        end
    end
end
