using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "Orthotope" begin
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

@testset "Simplex" begin
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

@testset "Subdivision" begin
    @testset "segment" begin
        s = HAI.segment(-0.81, 0.48)
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_segment2))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_segment3))
    end

    @testset "triangle" begin
        s = HAI.triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_triangle2))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_triangle4))
    end

    @testset "rectangle" begin
        s = HAI.rectangle((-0.07, 0.42), (0.35, 0.71))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_rectangle4))
    end

    @testset "tetrahedron" begin
        s = HAI.tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_tetrahedron8))
    end

    @testset "cuboid" begin
        s = HAI.cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test isnothing(HAI.check_subdivision(s, HAI.subdivide_cuboid8))
    end
end
