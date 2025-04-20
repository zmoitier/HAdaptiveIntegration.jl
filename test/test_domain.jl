using HAdaptiveIntegration:
    Cuboid,
    Orthotope,
    Rectangle,
    Segment,
    Simplex,
    Tetrahedron,
    Triangle,
    abs_det_jac,
    dimension,
    map_from_reference,
    map_to_reference,
    reference_domain,
    reference_orthotope,
    reference_simplex,
    subdivide_cuboid8,
    subdivide_orthotope,
    subdivide_rectangle4,
    subdivide_segment2,
    subdivide_simplex,
    subdivide_tetrahedron8,
    subdivide_triangle2,
    subdivide_triangle4,
    validate_subdivision
using StaticArrays
using Test

@testset "Domain construction" begin
    @testset "Orthotope" begin
        @test typeof(orthotope((0,), (1,))) <: Orthotope{1,Float64}
        @test typeof(orthotope([0], [1])) <: Orthotope{1,Float64}
        @test typeof(orthotope(SVector(0), SVector(1))) <: Orthotope{1,Float64}
        @test typeof(orthotope((big"0",), (big"1",))) <: Orthotope{1,BigFloat}

        r = reference_orthotope(Int, 4)
        @test typeof(r) <: Orthotope{4,Int}
        @test r.low_corner == SVector(0, 0, 0, 0)
        @test r.high_corner == SVector(1, 1, 1, 1)

        @test typeof(segment(-1, 1)) <: Orthotope{1,Float64}
        @test typeof(rectangle((-1, -1), (1, 1))) <: Orthotope{2,Float64}
        @test typeof(cuboid((-1, -1, -1), (1, 1, 1))) <: Orthotope{3,Float64}
        @test typeof(orthotope((-1, -1, -1, -1), (1, 1, 1, 1))) <: Orthotope{4,Float64}
    end

    @testset "Simplex" begin
        @test typeof(simplex((0,), (1,))) <: Simplex{1,2,Float64}
        @test typeof(simplex([0], [1])) <: Simplex{1,2,Float64}
        @test typeof(simplex(SVector(0), SVector(1))) <: Simplex{1,2,Float64}
        @test typeof(simplex((big"0",), (big"1",))) <: Simplex{1,2,BigFloat}

        r = reference_simplex(Int, 4)
        @test typeof(r) <: Simplex{4,5,Int}
        @test r.vertices[1] == SVector(0, 0, 0, 0)
        @test r.vertices[2] == SVector(1, 0, 0, 0)
        @test r.vertices[3] == SVector(0, 1, 0, 0)
        @test r.vertices[4] == SVector(0, 0, 1, 0)
        @test r.vertices[5] == SVector(0, 0, 0, 1)

        @test typeof(triangle((2, 0), (0, 2), (0, 0))) <: Simplex{2,3,Float64}
        @test typeof(tetrahedron((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))) <:
            Simplex{3,4,Float64}
        @test typeof(
            simplex((2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0))
        ) <: Simplex{4,5,Float64}
    end
end

@testset "Domain interface" begin
    @testset "Orthotope" begin
        h = orthotope((-1, -1, -1, -1), (1, 1, 1, 1))

        Φ = map_from_reference(h)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(-1, -1, -1, -1)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(-1, 1, -1, -1)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(-1, -1, 1, -1)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(-1, -1, -1, 1)

        Ψ = map_to_reference(h)
        @test Ψ(SVector(-1, -1, -1, -1)) ≈ SVector(0, 0, 0, 0)
        @test Ψ(SVector(-1, 1, -1, -1)) ≈ SVector(0, 1, 0, 0)
        @test Ψ(SVector(-1, -1, 1, -1)) ≈ SVector(0, 0, 1, 0)
        @test Ψ(SVector(-1, -1, -1, 1)) ≈ SVector(0, 0, 0, 1)

        @test abs_det_jac(h) ≈ 16

        @test typeof(reference_domain(Segment)) <: Orthotope{1,Float64}
        @test typeof(reference_domain(Segment{Int})) <: Orthotope{1,Int}
        @test typeof(reference_domain(Rectangle)) <: Orthotope{2,Float64}
        @test typeof(reference_domain(Rectangle{Int})) <: Orthotope{2,Int}
        @test typeof(reference_domain(Cuboid)) <: Orthotope{3,Float64}
        @test typeof(reference_domain(Cuboid{Int})) <: Orthotope{3,Int}
        @test typeof(reference_domain(Orthotope{4})) <: Orthotope{4,Float64}
        @test typeof(reference_domain(Orthotope{4,Int})) <: Orthotope{4,Int}

        @test dimension(Rectangle) == 2
        @test dimension(Rectangle{Int}) == 2
        @test dimension(Orthotope{4}) == 4
        @test dimension(Orthotope{4,Int}) == 4
    end

    @testset "Simplex" begin
        s = simplex((2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0))

        Φ = map_from_reference(s)
        @test Φ(SVector(0, 0, 0, 0)) ≈ SVector(2, 0, 0, 0)
        @test Φ(SVector(1, 0, 0, 0)) ≈ SVector(0, 2, 0, 0)
        @test Φ(SVector(0, 1, 0, 0)) ≈ SVector(0, 0, 2, 0)
        @test Φ(SVector(0, 0, 1, 0)) ≈ SVector(0, 0, 0, 2)
        @test Φ(SVector(0, 0, 0, 1)) ≈ SVector(0, 0, 0, 0)

        Ψ = map_to_reference(s)
        @test Ψ(SVector(2, 0, 0, 0)) ≈ SVector(0, 0, 0, 0)
        @test Ψ(SVector(0, 2, 0, 0)) ≈ SVector(1, 0, 0, 0)
        @test Ψ(SVector(0, 0, 2, 0)) ≈ SVector(0, 1, 0, 0)
        @test Ψ(SVector(0, 0, 0, 2)) ≈ SVector(0, 0, 1, 0)
        @test Ψ(SVector(0, 0, 0, 0)) ≈ SVector(0, 0, 0, 1)

        @test abs_det_jac(s) ≈ 16

        @test typeof(reference_domain(Triangle)) <: Simplex{2,3,Float64}
        @test typeof(reference_domain(Triangle{Int})) <: Simplex{2,3,Int}
        @test typeof(reference_domain(Tetrahedron)) <: Simplex{3,4,Float64}
        @test typeof(reference_domain(Tetrahedron{Int})) <: Simplex{3,4,Int}
        @test typeof(reference_domain(Simplex{4})) <: Simplex{4,5,Float64}
        @test typeof(reference_domain(Simplex{4,5})) <: Simplex{4,5,Float64}
        @test typeof(reference_domain(Simplex{4,5,Float64})) <: Simplex{4,5,Float64}

        @test dimension(Triangle) == 2
        @test dimension(Triangle{Int}) == 2
        @test dimension(Simplex{4}) == 4
        @test dimension(Simplex{4,5}) == 4
        @test dimension(Simplex{4,5,Int}) == 4
    end
end

@testset "Subdivision" begin
    @testset "segment" begin
        s = segment(-0.81, 0.48)
        @test validate_subdivision(subdivide_segment2, s)
    end

    @testset "triangle" begin
        t = triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test validate_subdivision(subdivide_triangle2, t)
        @test validate_subdivision(subdivide_triangle4, t)
        @test validate_subdivision(subdivide_simplex, t)
    end

    @testset "rectangle" begin
        r = rectangle((-0.07, 0.42), (0.35, 0.71))
        @test validate_subdivision(subdivide_rectangle4, r)
        @test validate_subdivision(subdivide_orthotope, r)
    end

    @testset "tetrahedron" begin
        t = tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test validate_subdivision(subdivide_tetrahedron8, t)
        @test validate_subdivision(subdivide_simplex, t)
    end

    @testset "cuboid" begin
        c = cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test validate_subdivision(subdivide_cuboid8, c)
        @test validate_subdivision(subdivide_orthotope, c)
    end

    @testset "4-simplex" begin
        s = simplex(
            (-0.13, -0.78, -0.22, 0.0),
            (0.70, -0.23, -0.37, 0.0),
            (-0.06, 0.57, -0.34, 0.0),
            (-0.27, -0.12, 0.14, 0.0),
            (-0.27, -0.12, 0.14, 1.0),
        )
        @test validate_subdivision(subdivide_simplex, s)
    end

    @testset "4-orthotope" begin
        h = orthotope((-2.12, -0.37, -0.86, 0.09), (-1.57, 0.11, 0.49, 0.66))
        @test validate_subdivision(subdivide_orthotope, h)
    end
end
