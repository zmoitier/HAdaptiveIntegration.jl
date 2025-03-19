using Test, StaticArrays
import HAdaptiveIntegration as hai

@testset "Domain construction" begin
    @testset "Orthotope" begin
        @test typeof(hai.orthotope((0,), (1,))) <: hai.Orthotope{1,Float64}
        @test typeof(hai.orthotope([0], [1])) <: hai.Orthotope{1,Float64}
        @test typeof(hai.orthotope(SVector(0), SVector(1))) <: hai.Orthotope{1,Float64}
        @test typeof(hai.orthotope((big"0",), (big"1",))) <: hai.Orthotope{1,BigFloat}

        r = hai.reference_orthotope(Int, 4)
        @test typeof(r) <: hai.Orthotope{4,Int}
        @test r.low_corner == SVector(0, 0, 0, 0)
        @test r.high_corner == SVector(1, 1, 1, 1)

        @test typeof(hai.segment(-1, 1)) <: hai.Orthotope{1,Float64}
        @test typeof(hai.rectangle((-1, -1), (1, 1))) <: hai.Orthotope{2,Float64}
        @test typeof(hai.cuboid((-1, -1, -1), (1, 1, 1))) <: hai.Orthotope{3,Float64}
        @test typeof(hai.orthotope((-1, -1, -1, -1), (1, 1, 1, 1))) <:
            hai.Orthotope{4,Float64}
    end

    @testset "Simplex" begin
        @test typeof(hai.simplex((0,), (1,))) <: hai.Simplex{1,2,Float64}
        @test typeof(hai.simplex([0], [1])) <: hai.Simplex{1,2,Float64}
        @test typeof(hai.simplex(SVector(0), SVector(1))) <: hai.Simplex{1,2,Float64}
        @test typeof(hai.simplex((big"0",), (big"1",))) <: hai.Simplex{1,2,BigFloat}

        r = hai.reference_simplex(Int, 4)
        @test typeof(r) <: hai.Simplex{4,5,Int}
        @test r.vertices[1] == SVector(0, 0, 0, 0)
        @test r.vertices[2] == SVector(1, 0, 0, 0)
        @test r.vertices[3] == SVector(0, 1, 0, 0)
        @test r.vertices[4] == SVector(0, 0, 1, 0)
        @test r.vertices[5] == SVector(0, 0, 0, 1)

        @test typeof(hai.triangle((2, 0), (0, 2), (0, 0))) <: hai.Simplex{2,3,Float64}
        @test typeof(hai.tetrahedron((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))) <:
            hai.Simplex{3,4,Float64}
        @test typeof(
            hai.simplex(
                (2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0)
            ),
        ) <: hai.Simplex{4,5,Float64}
    end
end

@testset "Domain interface" begin
    @testset "Orthotope" begin
        h = hai.orthotope((-1, -1, -1, -1), (1, 1, 1, 1))

        Φ = hai.map_from_reference(h)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(-1, -1, -1, -1)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(-1, 1, -1, -1)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(-1, -1, 1, -1)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(-1, -1, -1, 1)

        Ψ = hai.map_to_reference(h)
        @test Ψ(SVector(-1, -1, -1, -1)) ≈ SVector(0, 0, 0, 0)
        @test Ψ(SVector(-1, 1, -1, -1)) ≈ SVector(0, 1, 0, 0)
        @test Ψ(SVector(-1, -1, 1, -1)) ≈ SVector(0, 0, 1, 0)
        @test Ψ(SVector(-1, -1, -1, 1)) ≈ SVector(0, 0, 0, 1)

        @test hai.abs_det_jac(h) ≈ 16

        @test typeof(hai.reference_domain(hai.Segment)) <: hai.Orthotope{1,Float64}
        @test typeof(hai.reference_domain(hai.Segment{Int})) <: hai.Orthotope{1,Int}
        @test typeof(hai.reference_domain(hai.Rectangle)) <: hai.Orthotope{2,Float64}
        @test typeof(hai.reference_domain(hai.Rectangle{Int})) <: hai.Orthotope{2,Int}
        @test typeof(hai.reference_domain(hai.Cuboid)) <: hai.Orthotope{3,Float64}
        @test typeof(hai.reference_domain(hai.Cuboid{Int})) <: hai.Orthotope{3,Int}
        @test typeof(hai.reference_domain(hai.Orthotope{4})) <: hai.Orthotope{4,Float64}
        @test typeof(hai.reference_domain(hai.Orthotope{4,Int})) <: hai.Orthotope{4,Int}
    end

    @testset "Simplex" begin
        s = hai.simplex(
            (2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0)
        )

        Φ = hai.map_from_reference(s)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(2, 0, 0, 0)
        @test Φ(SVector(1, 0, 0, 0)) ≈ SVector(0, 2, 0, 0)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(0, 0, 2, 0)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(0, 0, 0, 2)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(0, 0, 0, 0)

        Ψ = hai.map_to_reference(s)
        @test Ψ(SVector(2, 0, 0, 0)) ≈ SVector(0, 0, 0, 0)
        @test Ψ(SVector(0, 2, 0, 0)) ≈ SVector(1, 0, 0, 0)
        @test Ψ(SVector(0, 0, 2, 0)) ≈ SVector(0, 1, 0, 0)
        @test Ψ(SVector(0, 0, 0, 2)) ≈ SVector(0, 0, 1, 0)
        @test Ψ(SVector(0, 0, 0, 0)) ≈ SVector(0, 0, 0, 1)

        @test hai.abs_det_jac(s) ≈ 16

        @test typeof(hai.reference_domain(hai.Triangle)) <: hai.Simplex{2,3,Float64}
        @test typeof(hai.reference_domain(hai.Triangle{Int})) <: hai.Simplex{2,3,Int}
        @test typeof(hai.reference_domain(hai.Tetrahedron)) <: hai.Simplex{3,4,Float64}
        @test typeof(hai.reference_domain(hai.Tetrahedron{Int})) <: hai.Simplex{3,4,Int}
        @test typeof(hai.reference_domain(hai.Simplex{4})) <: hai.Simplex{4,5,Float64}
        @test typeof(hai.reference_domain(hai.Simplex{4,Int})) <: hai.Simplex{4,5,Int}
        @test typeof(hai.reference_domain(hai.Simplex{4,5,Float64})) <:
            hai.Simplex{4,5,Float64}
    end
end

@testset "Subdivision" begin
    @testset "segment" begin
        s = hai.segment(-0.81, 0.48)
        @test hai.check_subdivision(hai.subdivide_segment2, s) == 0
    end

    @testset "triangle" begin
        t = hai.triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test hai.check_subdivision(hai.subdivide_triangle2, t) == 0
        @test hai.check_subdivision(hai.subdivide_triangle4, t) == 0
        @test hai.check_subdivision(hai.subdivide_simplex, t) == 0
    end

    @testset "rectangle" begin
        r = hai.rectangle((-0.07, 0.42), (0.35, 0.71))
        @test hai.check_subdivision(hai.subdivide_rectangle4, r) == 0
    end

    @testset "tetrahedron" begin
        t = hai.tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test hai.check_subdivision(hai.subdivide_tetrahedron8, t) == 0
        @test hai.check_subdivision(hai.subdivide_simplex, t) == 0
    end

    @testset "4-simplex" begin
        s = hai.simplex(
            (-0.13, -0.78, -0.22, 0.0),
            (0.70, -0.23, -0.37, 0.0),
            (-0.06, 0.57, -0.34, 0.0),
            (-0.27, -0.12, 0.14, 0.0),
            (-0.27, -0.12, 0.14, 1.0),
        )
        @test hai.check_subdivision(hai.subdivide_simplex, s) == 0
    end

    @testset "cuboid" begin
        c = hai.cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test hai.check_subdivision(hai.subdivide_cuboid8, c) == 0
    end
end
