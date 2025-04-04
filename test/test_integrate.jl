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
    segment = hai.segment(0, 1)
    buffer = hai.allocate_buffer(x -> zero(x[1]), segment)

    for (fct, R) in [
        (x -> exp(x[1]), exp(1) - 1),
        (x -> cos(10 * x[1]), sin(10) / 10),
        (x -> 1 / √x[1], 2),
    ]
        I, E = hai.integrate(fct, segment; buffer=buffer)
        @test abs(I - R) ≤ E * abs(R)
    end
end

@testset "Integrate over a triangle" begin
    triangle = hai.triangle((0, 0), (2, 0), (0, 2))
    buffer = hai.allocate_buffer(x -> zero(x[1]), triangle)

    for ec in
        (hai.default_embedded_cubature(triangle), hai.embedded_cubature(hai.TRIANGLE_GM19))
        for (fct, R) in [
            (x -> exp(x[1] + 3 * x[2]), (exp(6) - 3 * exp(2) + 2) / 6),
            (x -> cos(7 * x[1] + 3 * x[2]), (-3 * cos(14) + 7 * cos(6) - 4) / 84),
            (x -> 1 / norm(x), 2 * sqrt(2) * asinh(1)),
        ]
            I, E = hai.integrate(fct, triangle; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a rectangle" begin
    square = hai.rectangle((0, 0), (1, 1))
    buffer = hai.allocate_buffer(x -> zero(x[1]), square)

    for ec in (
        hai.default_embedded_cubature(square),
        hai.embedded_cubature(hai.SQUARE_CH21),
        hai.embedded_cubature(hai.SQUARE_GM17),
    )
        for (fct, R) in [
            (x -> exp(x[1] + x[2]), (exp(1) - 1)^2),
            (x -> 1 / norm(x), log(17 + 12 * sqrt(2)) / 2),
        ]
            I, E = hai.integrate(fct, square; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a tetrahedron" begin
    triangle = hai.tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    buffer = hai.allocate_buffer(x -> zero(x[1]), triangle)

    for (fct, R) in [
        (
            x -> exp(x[1] + 3 * x[2] + 5 * x[3]),
            (3 * exp(5) - 10 * exp(3) + 15 * exp(1) - 8) / 120,
        ),
        (x -> 1 / norm(x), 0.3614258523411),
    ]
        I, E = hai.integrate(fct, triangle; buffer=buffer)
        @test abs(I - R) ≤ E * abs(R)
    end
end

@testset "Integrate over a Cuboid" begin
    cube = hai.cuboid((0, 0, 0), (1, 1, 1))
    buffer = hai.allocate_buffer(x -> zero(x[1]), cube)

    for ec in (hai.default_embedded_cubature(cube), hai.embedded_cubature(hai.CUBE_GM33))
        for (fct, R) in [
            (x -> exp(x[1]), exp(1) - 1),
            (x -> 1 / (1 + norm(x)^2)^2, π^2 / 32),
            (x -> 1 / norm(x), 1.1900386819897766),
        ]
            I, E = hai.integrate(fct, cube; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "4-Simplex" begin
    I, E = hai.integrate(x -> 1 / norm(x), hai.reference_simplex(4))
    R = 0.089876019011
    @test abs(I - R) ≤ E * abs(R)
end

@testset "4-Orthotope" begin
    I, E = hai.integrate(x -> 1 / norm(x), hai.reference_orthotope(4))
    R = 0.9674120212411487
    @test abs(I - R) ≤ E * abs(R)
end
