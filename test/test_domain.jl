using HAdaptiveIntegration.Domain
using Quadmath
using StaticArrays
using Test

# Volume test
function volume_check(subdiv_algo, domain::AbstractDomain{D,T}) where {D,T}
    sub_domains = subdiv_algo(domain)
    return isapprox(sum(abs_det_jac.(sub_domains)), abs_det_jac(domain); rtol=10 * eps(T))
end

@testset "Domain construction" begin
    @testset "Simplex" begin
        @test typeof(Simplex((0, 0), [0, 0], SVector(0, 0))) <: Simplex{2,Int,3}
        @test typeof(Simplex((0, 0.0), [0, 0], SVector(0, 0))) <: Simplex{2,Float64,3}

        @test typeof(Simplex{Float64}([0], [0])) <: Simplex{1,Float64,2}
        @test typeof(Simplex{Float128}([0], [0])) <: Simplex{1,Float128,2}

        @test typeof(Triangle((2, 0), (0, 2), (0, 0))) <: Triangle{Int}
        @test typeof(Triangle{Float128}((2, 0), (0, 2), (0, 0))) <: Triangle{Float128}

        @test typeof(Tetrahedron((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))) <:
            Tetrahedron{Int}
        @test typeof(Tetrahedron{Float128}((0, 0, 2), (2, 0, 0), (0, 2, 0), (0, 0, 0))) <:
            Tetrahedron{Float128}
    end

    @testset "Orthotope" begin
        @test typeof(Orthotope((0,), [1])) <: Orthotope{1,Int}
        @test typeof(Orthotope([0], SVector(1.0))) <: Orthotope{1,Float64}

        @test typeof(Orthotope{Float64}([0], [1])) <: Orthotope{1,Float64}
        @test typeof(Orthotope{Float128}([0], [1])) <: Orthotope{1,Float128}

        @test typeof(Rectangle((-1, -1), (1, 1))) <: Rectangle{Int}
        @test typeof(Rectangle{Float128}((-1, -1), (1, 1))) <: Rectangle{Float128}

        @test typeof(Cuboid((-1, -1, -1), (1, 1, 1))) <: Cuboid{Int}
        @test typeof(Cuboid{Float128}((-1, -1, -1), (1, 1, 1))) <: Cuboid{Float128}
    end
end

@testset "Domain interface" begin
    @testset "Simplex" begin
        @test typeof(reference_domain(Triangle)) <: Triangle{Float64}
        @test typeof(reference_domain(Triangle{Float128})) <: Triangle{Float128}

        @test typeof(reference_domain(Tetrahedron)) <: Tetrahedron{Float64}
        @test typeof(reference_domain(Tetrahedron{Float128})) <: Tetrahedron{Float128}

        @test typeof(reference_domain(Simplex{4})) <: Simplex{4,Float64,5}
        @test typeof(reference_domain(Simplex{4,Float128})) <: Simplex{4,Float128,5}
        @test typeof(reference_domain(Simplex{4,Float128,5})) <: Simplex{4,Float128,5}

        r = reference_domain(Simplex{4,Int})
        @test r.vertices == SVector(
            SVector(0, 0, 0, 0),
            SVector(1, 0, 0, 0),
            SVector(0, 1, 0, 0),
            SVector(0, 0, 1, 0),
            SVector(0, 0, 0, 1),
        )

        s = Simplex(
            (-0.13, -0.78, -0.22, 0.04),
            (0.70, -0.23, -0.37, -0.72),
            (-0.06, 0.57, -0.34, 1.05),
            (-0.27, -0.12, 0.14, -0.47),
            (0.68, 0.12, -0.66, -1.49),
        )
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
        @test dimension(Simplex{4,Int,5}) == 4
    end

    @testset "Orthotope" begin
        @test typeof(reference_domain(Rectangle)) <: Orthotope{2,Float64}
        @test typeof(reference_domain(Rectangle{Float128})) <: Orthotope{2,Float128}

        @test typeof(reference_domain(Cuboid)) <: Orthotope{3,Float64}
        @test typeof(reference_domain(Cuboid{Float128})) <: Orthotope{3,Float128}

        @test typeof(reference_domain(Orthotope{4})) <: Orthotope{4,Float64}
        @test typeof(reference_domain(Orthotope{4,Float128})) <: Orthotope{4,Float128}

        r = reference_domain(Orthotope{4,Int})
        @test typeof(r) <: Orthotope{4,Int}
        @test r.corners == SVector(SVector(0, 0, 0, 0), SVector(1, 1, 1, 1))

        h = Orthotope((0.21, -0.58, -0.98, -1.25), (0.43, 1.75, 0.65, 1.87))
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
