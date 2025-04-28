using HAdaptiveIntegration.Domain
using StaticArrays
using Test

# Volume test
function volume_check(subdiv_algo, domain::AbstractDomain{D,T}) where {D,T}
    sub_domains = subdiv_algo(domain)
    return isapprox(sum(abs_det_jac.(sub_domains)), abs_det_jac(domain); rtol=10 * eps(T))
end

@testset "Domain construction" begin
    @testset "Simplex" begin
        @test typeof(Simplex((0,), (1,))) <: Simplex{1,Float64,2}
        @test typeof(Simplex([0], [1])) <: Simplex{1,Float64,2}
        @test typeof(Simplex(SVector(0.0), SVector(1.0))) <: Simplex{1,Float64,2}
        @test typeof(Simplex((big"0",), (big"1",))) <: Simplex{1,BigFloat,2}

        r = reference_domain(Simplex{4,Int})
        @test typeof(r) <: Simplex{4,Int,5}
        @test r.vertices[1] == SVector(0, 0, 0, 0)
        @test r.vertices[2] == SVector(1, 0, 0, 0)
        @test r.vertices[3] == SVector(0, 1, 0, 0)
        @test r.vertices[4] == SVector(0, 0, 1, 0)
        @test r.vertices[5] == SVector(0, 0, 0, 1)

        @test typeof(Triangle((2, 0), (0, 2), (0, 0))) <: Simplex{2,Float64,3}
        @test typeof(Tetrahedron((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))) <:
            Simplex{3,Float64,4}
        @test typeof(
            Simplex((2, 0, 0, 0), (0, 2, 0, 0), (0, 0, 2, 0), (0, 0, 0, 2), (0, 0, 0, 0))
        ) <: Simplex{4,Float64,5}
    end

    @testset "Orthotope" begin
        @test typeof(Orthotope((0,), (1,))) <: Orthotope{1,Float64}
        @test typeof(Orthotope([0], [1])) <: Orthotope{1,Float64}
        @test typeof(Orthotope(SVector(0.0), SVector(1.0))) <: Orthotope{1,Float64}
        @test typeof(Orthotope((big"0",), (big"1",))) <: Orthotope{1,BigFloat}

        r = reference_domain(Orthotope{4,Int})
        @test typeof(r) <: Orthotope{4,Int}
        @test r.low_corner == SVector(0, 0, 0, 0)
        @test r.high_corner == SVector(1, 1, 1, 1)

        @test typeof(Rectangle((-1, -1), (1, 1))) <: Orthotope{2,Float64}
        @test typeof(Cuboid((-1, -1, -1), (1, 1, 1))) <: Orthotope{3,Float64}
        @test typeof(Orthotope((-1, -1, -1, -1), (1, 1, 1, 1))) <: Orthotope{4,Float64}
    end
end

@testset "Domain interface" begin
    @testset "Simplex" begin
        s = Simplex(
            (-0.13, -0.78, -0.22, 0.04),
            (0.70, -0.23, -0.37, -0.72),
            (-0.06, 0.57, -0.34, 1.05),
            (-0.27, -0.12, 0.14, -0.47),
            (0.68, 0.12, -0.66, -1.49),
        )

        @test typeof(reference_domain(Triangle)) <: Simplex{2,Float64,3}
        @test typeof(reference_domain(Triangle{Int})) <: Simplex{2,Int,3}
        @test typeof(reference_domain(Tetrahedron)) <: Simplex{3,Float64,4}
        @test typeof(reference_domain(Tetrahedron{Int})) <: Simplex{3,Int,4}
        @test typeof(reference_domain(Simplex{4})) <: Simplex{4,Float64,5}
        @test typeof(reference_domain(Simplex{4,Int})) <: Simplex{4,Int,5}
        @test typeof(reference_domain(Simplex{4,Int,5})) <: Simplex{4,Int,5}

        Φ = map_from_reference(s)
        Ψ = map_to_reference(s)
        for v in (
            SVector(0, 0, 0, 0),
            SVector(1, 0, 0, 0),
            SVector(0, 1, 0, 0),
            SVector(0, 0, 1, 0),
            SVector(0, 0, 0, 1),
        )
            @test Ψ(Φ(v)) ≈ v
        end

        @test abs_det_jac(reference_domain(Simplex{4,Int})) ≈ 1

        @test dimension(Triangle) == 2
        @test dimension(Tetrahedron) == 3
        @test dimension(Simplex{4}) == 4
        @test dimension(Simplex{4,Int}) == 4
    end

    @testset "Orthotope" begin
        h = Orthotope((0.21, -0.58, -0.98, -1.25), (0.43, 1.75, 0.65, 1.87))

        @test typeof(reference_domain(Rectangle)) <: Orthotope{2,Float64}
        @test typeof(reference_domain(Rectangle{Int})) <: Orthotope{2,Int}
        @test typeof(reference_domain(Cuboid)) <: Orthotope{3,Float64}
        @test typeof(reference_domain(Cuboid{Int})) <: Orthotope{3,Int}
        @test typeof(reference_domain(Orthotope{4})) <: Orthotope{4,Float64}
        @test typeof(reference_domain(Orthotope{4,Int})) <: Orthotope{4,Int}

        Φ = map_from_reference(h)
        Ψ = map_to_reference(h)
        for v in (
            SVector(0, 0, 0, 0),
            SVector(1, 0, 0, 0),
            SVector(0, 1, 0, 0),
            SVector(0, 0, 1, 0),
            SVector(0, 0, 0, 1),
        )
            @test Ψ(Φ(v)) ≈ v
        end

        @test abs_det_jac(reference_domain(Orthotope{4,Int})) ≈ 1

        @test dimension(Rectangle) == 2
        @test dimension(Cuboid) == 3
        @test dimension(Orthotope{4}) == 4
        @test dimension(Orthotope{4,Int}) == 4
    end
end

@testset "Domain Subdivision" begin
    @testset "triangle" begin
        t = Triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test volume_check(subdivide_triangle, t)
        @test volume_check(subdivide_simplex, t)
    end

    @testset "rectangle" begin
        r = Rectangle((-0.07, 0.42), (0.35, 0.71))
        @test volume_check(subdivide_rectangle, r)
        @test volume_check(subdivide_orthotope, r)
    end

    @testset "tetrahedron" begin
        t = Tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test volume_check(subdivide_tetrahedron, t)
        @test volume_check(subdivide_simplex, t)
    end

    @testset "cuboid" begin
        c = Cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test volume_check(subdivide_cuboid, c)
        @test volume_check(subdivide_orthotope, c)
    end

    @testset "4-simplex" begin
        s = Simplex(
            (-0.13, -0.78, -0.22, 0.04),
            (0.70, -0.23, -0.37, -0.72),
            (-0.06, 0.57, -0.34, 1.05),
            (-0.27, -0.12, 0.14, -0.47),
            (0.68, 0.12, -0.66, -1.49),
        )
        @test volume_check(subdivide_simplex, s)
    end

    @testset "4-orthotope" begin
        h = Orthotope((-2.12, -0.37, -0.86, 0.09), (-1.57, 0.11, 0.49, 0.66))
        @test volume_check(subdivide_orthotope, h)
    end
end
