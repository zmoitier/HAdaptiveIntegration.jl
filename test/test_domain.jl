using HAdaptiveIntegration.Domain
using Quadmath
using StaticArrays
using Test

Float = float(Int)

# Volume test
function volume_test(subdiv_algo, domain::AbstractDomain{D,T}) where {D,T}
    sub_domains = subdiv_algo(domain)
    return isapprox(sum(abs_det_jac.(sub_domains)), abs_det_jac(domain); rtol=10 * eps(T))
end

@testset "Domain construction" begin
    @testset "Simplex" begin
        @test typeof(Segment(0, 1)) <: Segment{Float}
        @test typeof(Segment{Float128}(0, 1)) <: Segment{Float128}
        @test typeof(Segment(BigInt(0), BigInt(1))) <: Segment{BigFloat}
    end

    @testset "Simplex" begin
        @test typeof(Simplex((0, 0), [0, 0], SVector(0, 0))) <: Simplex{2,Float,3}
        @test typeof(Simplex{Float128}([0], [0])) <: Simplex{1,Float128,2}
        @test typeof(Simplex(BigInt[0], BigInt[0])) <: Simplex{1,BigFloat,2}

        @test typeof(Triangle((0, 0), (0, 0), (0, 0))) <: Triangle{Float}
        @test typeof(Triangle{Float128}((0, 0), (0, 0), (0, 0))) <: Triangle{Float128}
        @test typeof(Triangle(BigInt[0, 0], BigInt[0, 0], BigInt[0, 0])) <:
            Triangle{BigFloat}

        @test typeof(
            Tetrahedron(zeros(Int, 3), zeros(Int, 3), zeros(Int, 3), zeros(Int, 3))
        ) <: Tetrahedron{Float}
        @test typeof(
            Tetrahedron{Float128}(
                zeros(Int, 3), zeros(Int, 3), zeros(Int, 3), zeros(Int, 3)
            ),
        ) <: Tetrahedron{Float128}
        @test typeof(
            Tetrahedron(
                zeros(BigInt, 3), zeros(BigInt, 3), zeros(BigInt, 3), zeros(BigInt, 3)
            ),
        ) <: Tetrahedron{BigFloat}
    end

    @testset "Orthotope" begin
        @test typeof(Orthotope((0,), [1])) <: Orthotope{1,Float}
        @test typeof(Orthotope(BigInt[0], SVector(1))) <: Orthotope{1,BigFloat}
        @test typeof(Orthotope{Float128}([0], [1])) <: Orthotope{1,Float128}

        @test typeof(Rectangle((-1, -1), (1, 1))) <: Rectangle{Float}
        @test typeof(Rectangle{Float128}((-1, -1), (1, 1))) <: Rectangle{Float128}
        @test typeof(Rectangle(zeros(BigInt, 2), ones(BigInt, 2))) <: Rectangle{BigFloat}

        @test typeof(Cuboid((-1, -1, -1), (1, 1, 1))) <: Cuboid{Float}
        @test typeof(Cuboid{Float128}((-1, -1, -1), (1, 1, 1))) <: Cuboid{Float128}
        @test typeof(Cuboid(zeros(BigInt, 3), ones(BigInt, 3))) <: Cuboid{BigFloat}
    end
end

@testset "Domain interface" begin
    @testset "Segment" begin
        @test typeof(reference_domain(Segment)) <: Segment{Float}
        @test typeof(reference_domain(Segment{Float128})) <: Segment{Float128}
        @test typeof(reference_domain(Segment{BigFloat})) <: Segment{BigFloat}

        r = reference_domain(Segment{Int})
        @test r.xmin == 0
        @test r.xmax == 1

        s = Segment(-0.13, 0.78)
        Φ = map_from_reference(s)
        Ψ = map_to_reference(s)
        for v in (0, 1)
            @test Ψ(Φ(v)) ≈ v
        end

        @test abs_det_jac(r) ≈ 1

        @test dimension(Segment) == 1
        @test dimension(Segment{Int}) == 1
    end

    @testset "Simplex" begin
        @test typeof(reference_domain(Triangle)) <: Triangle{Float}
        @test typeof(reference_domain(Triangle{Float128})) <: Triangle{Float128}

        @test typeof(reference_domain(Tetrahedron)) <: Tetrahedron{Float}
        @test typeof(reference_domain(Tetrahedron{Float128})) <: Tetrahedron{Float128}

        @test typeof(reference_domain(Simplex{4})) <: Simplex{4,Float,5}
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
        @test typeof(reference_domain(Rectangle)) <: Orthotope{2,Float}
        @test typeof(reference_domain(Rectangle{Float128})) <: Orthotope{2,Float128}

        @test typeof(reference_domain(Cuboid)) <: Orthotope{3,Float}
        @test typeof(reference_domain(Cuboid{Float128})) <: Orthotope{3,Float128}

        @test typeof(reference_domain(Orthotope{4})) <: Orthotope{4,Float}
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
    @testset "segment" begin
        s = Segment(-0.13, 0.78)
        @test volume_test(subdivide_segment, s)
    end

    @testset "triangle" begin
        t = Triangle((-0.86, -0.19), (0.97, -0.84), (-0.05, 0.74))
        @test volume_test(subdivide_triangle, t)
        @test volume_test(subdivide_simplex, t)
    end

    @testset "rectangle" begin
        r = Rectangle((-0.07, 0.42), (0.35, 0.71))
        @test volume_test(subdivide_rectangle, r)
        @test volume_test(subdivide_orthotope, r)
    end

    @testset "tetrahedron" begin
        t = Tetrahedron(
            (-0.13, -0.78, -0.22),
            (0.70, -0.23, -0.37),
            (-0.06, 0.57, -0.34),
            (-0.27, -0.12, 0.14),
        )
        @test volume_test(subdivide_tetrahedron, t)
        @test volume_test(subdivide_simplex, t)
    end

    @testset "cuboid" begin
        c = Cuboid((-0.38, -0.92, -0.43), (0.15, 0.30, 0.51))
        @test volume_test(subdivide_cuboid, c)
        @test volume_test(subdivide_orthotope, c)
    end

    @testset "4-simplex" begin
        s = Simplex(
            (-0.13, -0.78, -0.22, 0.04),
            (0.70, -0.23, -0.37, -0.72),
            (-0.06, 0.57, -0.34, 1.05),
            (-0.27, -0.12, 0.14, -0.47),
            (0.68, 0.12, -0.66, -1.49),
        )
        @test volume_test(subdivide_simplex, s)
    end

    @testset "4-orthotope" begin
        h = Orthotope((-2.12, -0.37, -0.86, 0.09), (-1.57, 0.11, 0.49, 0.66))
        @test volume_test(subdivide_orthotope, h)
    end
end
