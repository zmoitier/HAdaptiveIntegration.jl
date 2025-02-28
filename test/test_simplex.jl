using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "domain_simplex.jl" begin
    @testset "Simplex construction" begin
        @test typeof(HAI.simplex((0,), (1,))) <: HAI.Simplex
        @test typeof(HAI.simplex([0], [1])) <: HAI.Simplex
        @test typeof(HAI.simplex(SVector{1}([0]), SVector{1}([1]))) <: HAI.Simplex
    end

    @testset "Triangle" begin
        ref = HAI.reference_domain(HAI.Triangle{Float64})
        t = HAI.triangle((0, 0), (1, 0), (0, 1))
        @test t.points == ref.points
        @test HAI.abs_det_jacobian(ref) ≈ 1
    end

    @testset "Tetrahedron" begin
        ref = HAI.reference_domain(HAI.Tetrahedron{Float64})
        t = HAI.tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        @test t.points == ref.points
        @test HAI.abs_det_jacobian(ref) ≈ 1
    end

    @testset "4-simplex" begin
        ref = HAI.reference_domain(HAI.Simplex{4,Float64,5})
        @test HAI.abs_det_jacobian(ref) ≈ 1

        Φ = HAI.map_from_reference(ref)
        @test Φ([0, 0, 0, 0]) ≈ SVector{4}([0, 0, 0, 0])
        @test Φ([0, 1, 0, 0]) ≈ SVector{4}([0, 1, 0, 0])
        @test Φ([0, 0, 1, 0]) ≈ SVector{4}([0, 0, 1, 0])
        @test Φ([0, 0, 0, 1]) ≈ SVector{4}([0, 0, 0, 1])
    end
end
