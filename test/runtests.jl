import AdaptiveSimplexQuadrature as ASQ
using Test
using StaticArrays
using LinearAlgebra

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
    @testset "Triangle" begin
        tri = ASQ.Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        @test ASQ.measure(tri) ≈ 0.5
        f = x -> 1
        quad = ASQ.EmbeddedQuadrature(; name = "LaurieRadon")
        I, E = ASQ._integrate_with_error(f, tri, quad)
        @test I ≈ 0.5
        counter = Int[0]
        f = x -> begin
            counter[1] += 1
            return log(norm(x))
        end
        I, E = ASQ.integrate(f, tri; atol = 1e-2)
        @show I, E
    end
end
