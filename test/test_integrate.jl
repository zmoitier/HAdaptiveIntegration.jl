using Test, StaticArrays, LinearAlgebra, DataStructures
import HAdaptiveIntegration as hai

@testset "Embedded cubature" begin
    tec = hai.TabulatedEmbeddedCubature(;
        description="Gauss (SEGMENT_G3)",
        domain="reference segment",
        reference="",
        nb_significant_digits=16,
        nodes=[["5e-1"], ["1.127016653792583e-1"], ["8.872983346207417e-1"]],
        weights_high=[
            "4.444444444444444e-1", "2.777777777777778e-1", "2.777777777777778e-1"
        ],
        order_high=5,
        weights_low=["1"],
        order_low=1,
    )
    @test typeof(tec) <: hai.TabulatedEmbeddedCubature

    ec = hai.embedded_cubature(tec, Float64)
    @test typeof(ec) <: hai.EmbeddedCubature{3,1,1,Float64}

    ec_ref = hai.embedded_cubature(
        [[0.5], [(1 - √(3 / 5)) / 2], [(1 + √(3 / 5)) / 2]], [4 / 9, 5 / 18, 5 / 18], [1.0]
    )
    @test typeof(ec_ref) <: hai.EmbeddedCubature{3,1,1,Float64}

    @test ec.nodes ≈ ec_ref.nodes
    @test ec.weights_high ≈ ec_ref.weights_high
    @test ec.weights_low ≈ ec_ref.weights_low

    @test typeof(hai.embedded_cubature(hai.GrundmannMoeller(2, 5))) <: hai.EmbeddedCubature
end

@testset "Default embedded cubature" begin
    for domain in (
        hai.reference_orthotope(1),
        hai.reference_orthotope(2),
        hai.reference_orthotope(3),
        hai.reference_simplex(2),
        hai.reference_simplex(3),
        hai.reference_simplex(4),
    )
        @test typeof(hai.default_embedded_cubature(domain)) <: hai.EmbeddedCubature
    end
end

@testset "Integrate" begin
    domain = hai.reference_orthotope(1)

    buffer = hai.allocate_buffer(x -> sum(x), domain)
    @test typeof(buffer) <: BinaryHeap

    I, E = hai.integrate(x -> sin(10 * x[1]), domain; buffer=buffer)
    R = sin(5)^2 / 5
    @test abs(I - R) ≤ E * abs(R)

    I, E = hai.integrate(x -> cos(7.5 * x[1]), domain; buffer=buffer)
    R = 2 * sin(7.5) / 15
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a segment" begin
    @test hai.check_order(hai.SEGMENT_GK15, hai.reference_orthotope(1)) == 0

    ec = hai.embedded_cubature(hai.SEGMENT_GK15, Float64)
    segment = hai.segment(0, 1)

    I, E = ec(x -> exp(x[1]), segment)
    R = exp(1) - exp(0)
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> cos(10 * x[1]), segment)
    R = sin(10) / 10
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> 1 / √x[1], segment)
    R = 2
    @test abs(I - R) ≤ E * abs(R)

    I, E = hai.integrate(x -> 1 / √x[1], segment)
    R = 2
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a rectangle" begin
    @test hai.check_order(hai.SQUARE_CHG25, hai.reference_orthotope(2)) == 0

    ec = hai.embedded_cubature(hai.SQUARE_CHG25, Float64)
    square = hai.rectangle((0, 0), (1, 1))

    I, E = ec(x -> exp(x[1] + x[2]), square)
    R = (exp(1) - exp(0))^2
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> 1 / norm(x), square)
    R = 1 / 2 * log(17 + 12 * sqrt(2))
    @test abs(I - R) ≤ E * abs(R)

    I, E = hai.integrate(x -> 1 / norm(x), square)
    @test abs(I - R) < 1e-8
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a triangle" begin
    @test hai.check_order(hai.TRIANGLE_RL19, hai.reference_simplex(2)) == 0

    ec = hai.embedded_cubature(hai.TRIANGLE_RL19, Float64)

    triangle = hai.triangle((0, 0), (2, 0), (0, 2))
    I, E = ec(x -> exp(x[1] + 3 * x[2]), triangle)
    R = (exp(6) - 3 * exp(2) + 2) / 6
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> cos(7 * x[1] + 3 * x[2]), triangle)
    R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
    @test abs(I - R) ≤ E * abs(R)

    triangle = hai.triangle((0, 0), (1, 0), (0, 1))
    I, E = hai.integrate(x -> 1 / norm(x), triangle)
    R = sqrt(2) * asinh(1)
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a Cuboid" begin
    @test hai.check_order(hai.CUBE_BE65, hai.reference_orthotope(3)) == 0

    ec = hai.embedded_cubature(hai.CUBE_BE65, Float64)
    cube = hai.cuboid((0, 0, 0), (1, 1, 1))

    I, E = ec(x -> exp(x[1]), cube)
    R = exp(1) - exp(0)
    @test abs(I - R) < 1e-8
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> exp(x[1] + x[2] + x[3]), cube)
    R = (exp(1) - exp(0))^3
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> 1 / (1 + norm(x)^2)^2, cube)
    R = π^2 / 32
    @test abs(I - R) ≤ E * abs(R)

    I, E = hai.integrate(x -> 1 / norm(x), cube)
    @test E ≤ √eps(Float64) * I
end

@testset "GrundmannMoeller quadrature" begin
    @testset "Triangle" begin
        @test hai.check_order(hai.TRIANGLE_GM20, hai.reference_simplex(2)) == 0

        ec = hai.embedded_cubature(hai.TRIANGLE_GM20, Float64)

        rtol = 1e-8

        triangle = hai.triangle((0, 0), (2, 0), (0, 2))
        I, E = hai.integrate(
            x -> exp(x[1] + 3 * x[2]), triangle; embedded_cubature=ec, rtol=rtol
        )
        R = (exp(6) - 3 * exp(2) + 2) / 6
        @test abs(I - R) ≤ rtol * abs(R)

        I, E = hai.integrate(
            x -> cos(7 * x[1] + 3 * x[2]), triangle; embedded_cubature=ec, rtol=rtol
        )
        R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
        @test abs(I - R) ≤ rtol * abs(R)

        triangle = hai.triangle((0, 0), (1, 0), (0, 1))
        I, E = hai.integrate(x -> 1 / norm(x), triangle; embedded_cubature=ec, rtol=rtol)
        R = sqrt(2) * asinh(1)
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "Tetrahedron" begin
        @test hai.check_order(
            hai.TETRAHEDRON_GM35, hai.reference_simplex(3); rtol=12 * eps(Float64)
        ) == 0

        tetrahedron = hai.tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        rtol = 1e-8

        I, E = hai.integrate(x -> 1, tetrahedron; rtol=rtol)
        R = 1 / 6
        @test abs(I - R) ≤ rtol * abs(R)

        I, E = hai.integrate(x -> exp(x[1] + 3 * x[2] + 5 * x[3]), tetrahedron; rtol=rtol)
        R = (3 * exp(5) - 10 * exp(3) + 15 * exp(1) - 8) / 120
        @test abs(I - R) ≤ rtol * abs(R)

        I, E = hai.integrate(x -> 1 / norm(x), tetrahedron; rtol=rtol)
        R = 0.3614258523411 # reference computed using BigFloat
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "4-Simplex" begin
        ec = hai.embedded_cubature(hai.GrundmannMoeller(4, 7), Float64)

        @test hai.check_order(ec, hai.reference_simplex(4), 7, 5; rtol=50 * eps(Float64)) ==
            0

        rtol = 1e-8
        I, E = hai.integrate(x -> 1 / norm(x), hai.reference_simplex(4); rtol=rtol)
        R = 0.089876019011 # reference computed using BigFloat
        @test abs(I - R) ≤ rtol * abs(R)
    end
end
