using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "domain_simplex.jl" begin
    @testset "Simplex construction" begin
        @test typeof(HAI.Simplex((0,), (1,))) <: HAI.Simplex
        @test typeof(HAI.Simplex([0], [1])) <: HAI.Simplex
        @test typeof(HAI.Simplex(SVector{1}([0]), SVector{1}([1]))) <: HAI.Simplex
    end

    @testset "Triangle" begin
        ref = HAI.reference_triangle()
        t = HAI.Triangle((0, 0), (1, 0), (0, 1))
        @test t.points == ref.points
        @test HAI.abs_jacobian_determinant(ref) ≈ 1
    end

    @testset "Tetrahedron" begin
        ref = HAI.reference_tetrahedron()
        t = HAI.Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        @test t.points == ref.points
        @test HAI.abs_jacobian_determinant(ref) ≈ 1
    end

    @testset "4-simplex" begin
        ref = HAI.reference_simplex(4)
        @test HAI.abs_jacobian_determinant(ref) ≈ 1

        Φ = HAI.map_from_reference(ref)
        @test Φ([0, 0, 0, 0]) ≈ SVector{4}([0, 0, 0, 0])
        @test Φ([0, 1, 0, 0]) ≈ SVector{4}([0, 1, 0, 0])
        @test Φ([0, 0, 1, 0]) ≈ SVector{4}([0, 0, 1, 0])
        @test Φ([0, 0, 0, 1]) ≈ SVector{4}([0, 0, 0, 1])
    end
end
