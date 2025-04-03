using DataStructures
using LinearAlgebra
using Test

import HAdaptiveIntegration as hai

@testset "Default embedded cubature" begin
    for domain in (
        hai.reference_orthotope(1),
        hai.reference_orthotope(2),
        hai.reference_orthotope(3),
        hai.reference_orthotope(4),
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
    ec = hai.embedded_cubature(Float64, hai.SEGMENT_GK15)
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
    ec = hai.embedded_cubature(Float64, hai.SQUARE_CH25)
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
    ec = hai.embedded_cubature(Float64, hai.TRIANGLE_RL19)

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
    ec = hai.embedded_cubature(Float64, hai.CUBE_BE65)
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
        ec = hai.embedded_cubature(Float64, hai.TRIANGLE_GM19)
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
        tetrahedron = hai.tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        ec = hai.embedded_cubature(Float64, hai.TETRAHEDRON_GM35)
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
        rtol = 1e-8
        I, E = hai.integrate(x -> 1 / norm(x), hai.reference_simplex(4); rtol=rtol)
        R = 0.089876019011 # reference computed using BigFloat
        @test abs(I - R) ≤ rtol * abs(R)
    end
end
