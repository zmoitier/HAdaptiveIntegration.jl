import AdaptiveSimplexQuadrature as ASQ
using Test

@testset "AdaptiveSimplexQuadrature.jl" begin
    # include("aqua_test.jl")

    @testset "Segment" begin
        segment = ASQ.Segment((0.0), (1.0))
        @test ASQ.det_jac(segment) ≈ 1

        @test ASQ.check_order(ASQ.SEGMENT_G7, 13)
        @test ASQ.check_order(ASQ.SEGMENT_K15, 23)
        embd_quad = ASQ.EmbeddedQuadrature(; name = "segment-G7K15")

        simplex = ASQ.Simplex((0,), (1,))
        I, E = embd_quad(x -> exp(x[1]), simplex)
        R = exp(1) - exp(0)
        @test abs(I - R) ≤ E * abs(R)

        simplex = ASQ.Simplex((0,), (1,))
        I, E = embd_quad(x -> cos(10 * x[1]), simplex)
        R = sin(10) / 10
        @test abs(I - R) ≤ E * abs(R)

        simplex = ASQ.Simplex((0,), (1,))
        I, E = embd_quad(x -> 1 / √x[1], simplex)
        R = 2
        @test abs(I - R) ≤ E * abs(R)

        simplex = ASQ.Simplex((0.0,), (1.0,))
        I, E = ASQ.integrate(x -> 1 / √x[1], simplex, embd_quad)
        R = 2
        @test abs(I - R) ≤ E * abs(R)
    end

    @testset "Triangle" begin
        triangle = ASQ.Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        @test ASQ.det_jac(triangle) ≈ 1

        @test ASQ.check_order(ASQ.TRIANGLE_R5N7, 5)
        @test ASQ.check_order(ASQ.TRIANGLE_L8N19, 8)

        embd_quad = ASQ.EmbeddedQuadrature(; name = "triangle-LaurieRadon")

        simplex = ASQ.Simplex((0, 0), (2, 0), (0, 2))
        I, E = embd_quad(x -> exp(x[1] + 3 * x[2]), simplex)
        R = (exp(6) - 3 * exp(2) + 2) / 6
        @test abs(I - R) ≤ E * abs(R)

        simplex = ASQ.Simplex((0, 0), (2, 0), (0, 2))
        I, E = embd_quad(x -> cos(7 * x[1] + 3 * x[2]), simplex)
        R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
        @test abs(I - R) ≤ E * abs(R)

        # simplex = ASQ.Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        # I, E = embd_quad(x -> 1 / √x[1], simplex)
        # R = 4 / 3
        # @test abs(I - R) ≤ E * abs(R)

        # simplex = ASQ.Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        # I, E = ASQ.integrate(x -> 1 / √x[1], simplex, embd_quad)
        # R = 4 / 3
        # @test abs(I - R) ≤ E * abs(R)
    end

    # @testset "Simplex" begin
    #     @test_throws AssertionError ASQ.Simplex(7
    #         (0.0, 0.0),
    #         (1.0, 0.0),
    #         (0.0, 1.0, 0.0, 0.0),
    #     )
    #     @test_throws AssertionError ASQ.Simplex((0.0, 0.0), (1.0, 0.0))
    #     tri = ASQ.Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    #     tris = ASQ.subdivide(tri)
    #     @test length(tris) == 4
    #     s = ASQ.Simplex((0.0, 0.0, 1.0), (1.0, 0.0, 0.5), (0.0, 1.0, 0.3), (-0.1, 0.2, 0.3))
    # end
end
