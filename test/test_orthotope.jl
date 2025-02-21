using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "domain_orthotope.jl" begin
    @testset "Orthotope construction" begin
        @test typeof(HAI.Orthotope((0,), (1,))) <: HAI.Orthotope
        @test typeof(HAI.Orthotope([0], [1])) <: HAI.Orthotope
        @test typeof(HAI.Orthotope(SVector{1}([0]), SVector{1}([1]))) <: HAI.Orthotope
    end

    @testset "Segment" begin
        ref = HAI.reference_segment()
        s = HAI.Segment((0,), (1,))
        @test (s.low_corner ≈ ref.low_corner) && (s.high_corner ≈ ref.high_corner)
        @test HAI.abs_jacobian_determinant(ref) ≈ 1
    end

    @testset "Rectangle" begin
        ref = HAI.reference_rectangle()
        r = HAI.Rectangle((0, 0), (1, 1))
        @test (r.low_corner ≈ ref.low_corner) && (r.high_corner ≈ ref.high_corner)
        @test HAI.abs_jacobian_determinant(ref) ≈ 1
    end

    @testset "Cuboid" begin
        ref = HAI.reference_cuboid()
        c = HAI.Cuboid((0, 0, 0), (1, 1, 1))
        @test (c.low_corner ≈ ref.low_corner) && (c.high_corner ≈ ref.high_corner)
        @test HAI.abs_jacobian_determinant(ref) ≈ 1
    end

    @testset "4-orthotope" begin
        ref = HAI.reference_orthotope(4)
        @test HAI.abs_jacobian_determinant(ref) ≈ 1

        Φ = HAI.map_from_reference(ref)
        @test Φ([0, 0, 0, 0]) ≈ SVector{4}([0, 0, 0, 0])
        @test Φ([0, 1, 0, 0]) ≈ SVector{4}([0, 1, 0, 0])
        @test Φ([0, 0, 1, 0]) ≈ SVector{4}([0, 0, 1, 0])
        @test Φ([0, 0, 0, 1]) ≈ SVector{4}([0, 0, 0, 1])
    end
end
