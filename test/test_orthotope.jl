using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "domain_orthotope.jl" begin
    @testset "Orthotope construction" begin
        @test typeof(HAI.orthotope((0,), (1,))) <: HAI.Orthotope
        @test typeof(HAI.orthotope([0], [1])) <: HAI.Orthotope
        @test typeof(HAI.orthotope(SVector{1}([0]), SVector{1}([1]))) <: HAI.Orthotope
    end

    @testset "Segment" begin
        ref = HAI.reference_domain(HAI.Segment{Float64})
        s = HAI.segment(0, 1)
        @test (s.low_corner ≈ ref.low_corner) && (s.high_corner ≈ ref.high_corner)
        @test HAI.abs_det_jacobian(ref) ≈ 1
    end

    @testset "Rectangle" begin
        ref = HAI.reference_domain(HAI.Rectangle{Float64})
        r = HAI.rectangle((0, 0), (1, 1))
        @test (r.low_corner ≈ ref.low_corner) && (r.high_corner ≈ ref.high_corner)
        @test HAI.abs_det_jacobian(ref) ≈ 1
    end

    @testset "Cuboid" begin
        ref = HAI.reference_domain(HAI.Cuboid{Float64})
        c = HAI.cuboid((0, 0, 0), (1, 1, 1))
        @test (c.low_corner ≈ ref.low_corner) && (c.high_corner ≈ ref.high_corner)
        @test HAI.abs_det_jacobian(ref) ≈ 1
    end

    @testset "4-orthotope" begin
        ref = HAI.reference_domain(HAI.Orthotope{4,Float64})
        @test HAI.abs_det_jacobian(ref) ≈ 1

        Φ = HAI.map_from_reference(ref)
        @test Φ([0, 0, 0, 0]) ≈ SVector{4}([0, 0, 0, 0])
        @test Φ([0, 1, 0, 0]) ≈ SVector{4}([0, 1, 0, 0])
        @test Φ([0, 0, 1, 0]) ≈ SVector{4}([0, 0, 1, 0])
        @test Φ([0, 0, 0, 1]) ≈ SVector{4}([0, 0, 0, 1])
    end
end
