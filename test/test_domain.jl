using Test, StaticArrays
import HAdaptiveIntegration as HAI

@testset "Orthotope" begin
    @testset "Orthotope construction" begin
        @test typeof(HAI.orthotope((0,), (1,))) <: HAI.Orthotope
        @test typeof(HAI.orthotope([0], [1.0])) <: HAI.Orthotope
        @test typeof(HAI.orthotope(SVector(0), SVector(1))) <: HAI.Orthotope
    end

    @testset "Segment" begin
        s = HAI.segment(-1, 1)
        @test HAI.abs_det_jac(s) ≈ 2
    end

    @testset "Rectangle" begin
        r = HAI.rectangle((-1, -1), (1, 1))
        @test HAI.abs_det_jac(r) ≈ 4
    end

    @testset "Cuboid" begin
        c = HAI.cuboid((-1, -1, -1), (1, 1, 1))
        @test HAI.abs_det_jac(c) ≈ 8
    end

    @testset "4-orthotope" begin
        h = HAI.orthotope((-1, -1, -1, -1), (1, 1, 1, 1))
        @test HAI.abs_det_jac(h) ≈ 16

        Φ = HAI.map_from_reference(h)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(-1, -1, -1, -1)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(-1, 1, -1, -1)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(-1, -1, 1, -1)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(-1, -1, -1, 1)
    end
end

@testset "Simplex" begin
    @testset "Simplex construction" begin
        @test typeof(HAI.simplex((0,), (1,))) <: HAI.Simplex
        @test typeof(HAI.simplex([0], [1])) <: HAI.Simplex
        @test typeof(HAI.simplex(SVector(0), SVector(1))) <: HAI.Simplex
    end

    @testset "Triangle" begin
        t = HAI.triangle((2, 0), (0, 2), (0, 0))
        @test HAI.abs_det_jac(t) ≈ 4
    end

    @testset "Tetrahedron" begin
        t = HAI.tetrahedron((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))
        @test HAI.abs_det_jac(t) ≈ 8
    end

    @testset "4-simplex" begin
        s = HAI.simplex(
            (2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0)
        )
        @test HAI.abs_det_jac(s) ≈ 16

        Φ = HAI.map_from_reference(s)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(2, 0, 0, 0)
        @test Φ(SVector(1, 0, 0, 0)) ≈ SVector(0, 2, 0, 0)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(0, 0, 2, 0)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(0, 0, 0, 2)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(0, 0, 0, 0)
    end
end

@testset "Subdivision" begin
    @testset "segment" begin
        s = HAI.segment(-0.81, 0.48)
        @test HAI.check_subdivision(HAI.subdivide_segment2, s) == 0
    end

    @testset "triangle" begin
        t = HAI.triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test HAI.check_subdivision(HAI.subdivide_triangle2, t) == 0
        @test HAI.check_subdivision(HAI.subdivide_triangle4, t) == 0
        @test HAI.check_subdivision(HAI.subdivide_simplex, t) == 0
    end

    @testset "rectangle" begin
        r = HAI.rectangle((-0.07, 0.42), (0.35, 0.71))
        @test HAI.check_subdivision(HAI.subdivide_rectangle4, r) == 0
    end

    @testset "tetrahedron" begin
        t = HAI.tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test HAI.check_subdivision(HAI.subdivide_tetrahedron8, t) == 0
        @test HAI.check_subdivision(HAI.subdivide_simplex, t) == 0
    end

    @testset "4-simplex" begin
        s = HAI.simplex(
            (-0.13, -0.78, -0.22, 0.0),
            (0.70, -0.23, -0.37, 0.0),
            (-0.06, 0.57, -0.34, 0.0),
            (-0.27, -0.12, 0.14, 0.0),
            (-0.27, -0.12, 0.14, 1.0),
        )
        @test HAI.check_subdivision(HAI.subdivide_simplex, s) == 0
    end

    @testset "cuboid" begin
        c = HAI.cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test HAI.check_subdivision(HAI.subdivide_cuboid8, c) == 0
    end
end
