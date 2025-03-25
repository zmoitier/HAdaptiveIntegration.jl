using Test
import HAdaptiveIntegration as hai

@testset "Embedded cubature" begin
    tec = hai.TabulatedEmbeddedCubature{hai.Segment}(;
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
    @test typeof(tec) <: hai.TabulatedEmbeddedCubature{hai.Segment}

    ec = hai.embedded_cubature(tec)
    @test typeof(ec) <: hai.EmbeddedCubature{1,Float64}

    ec_ref = hai.embedded_cubature(
        [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]], [4 / 9, 5 / 18, 5 / 18], [1.0]
    )
    @test typeof(ec_ref) <: hai.EmbeddedCubature{1,Float64}

    @test ec.nodes ≈ ec_ref.nodes
    @test ec.weights_high ≈ ec_ref.weights_high
    @test ec.weights_low ≈ ec_ref.weights_low

    @test typeof(hai.embedded_cubature(hai.GrundmannMoeller{2}(5))) <:
        hai.EmbeddedCubature{2,Float64}
end

@testset "Tabulated rules" begin
    @testset "Segment" begin
        @test hai.check_order(hai.SEGMENT_GK7, hai.reference_orthotope(1)) == 0
        @test hai.check_order(hai.SEGMENT_GK15, hai.reference_orthotope(1)) == 0
        @test hai.check_order(hai.SEGMENT_GK31, hai.reference_orthotope(1)) == 0
    end

    @testset "Square" begin
        @test hai.check_order(hai.SQUARE_CHG25, hai.reference_orthotope(2)) == 0
        @test hai.check_order(hai.SQUARE_CHG21, hai.reference_orthotope(2)) == 0
    end

    @testset "Triangle" begin
        @test hai.check_order(hai.TRIANGLE_RL19, hai.reference_simplex(2)) == 0
        @test hai.check_order(hai.TRIANGLE_GM20, hai.reference_simplex(2)) == 0
    end

    @testset "Cube" begin
        @test hai.check_order(hai.CUBE_BE65, hai.reference_orthotope(3)) == 0
    end

    @testset "Tetrahedron" begin
        @test hai.check_order(
            hai.TETRAHEDRON_GM35, hai.reference_simplex(3); rtol=12 * eps(Float64)
        ) == 0
    end
end

@testset "Grundmann Möller" begin
    @test hai.check_order(
        hai.embedded_cubature(Float64, hai.GrundmannMoeller{4}(7)),
        hai.reference_simplex(4),
        7,
        5;
        rtol=50 * eps(Float64),
    ) == 0
end
