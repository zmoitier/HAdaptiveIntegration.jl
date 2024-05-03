import AdaptiveSimplexQuadrature as ASQ
using Test
using StaticArrays

@testset "AdaptiveSimplexQuadrature.jl" begin
    include("aqua_test.jl")

    @testset "Simplex" begin
        @test_throws AssertionError ASQ.Simplex(
            (0.0, 0.0),
            (1.0, 0.0),
            (0.0, 1.0, 0.0, 0.0),
        )
        @test_throws AssertionError ASQ.Simplex((0.0, 0.0), (1.0, 0.0))
        tri = ASQ.Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        tris = ASQ.subdivide(tri)
        @test length(tris) == 4
        s = ASQ.Simplex((0.0, 0.0, 1.0), (1.0, 0.0, 0.5), (0.0, 1.0, 0.3), (-0.1, 0.2, 0.3))
    end
end
