using DataStructures: BinaryHeap
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using HAdaptiveIntegration: allocate_buffer, default_embedded_cubature, integrate
using LinearAlgebra: norm
using Test

@testset "Default embedded cubature" begin
    for domain in
        (Segment, Rectangle, Cuboid, Orthotope{4}, Triangle, Tetrahedron, Simplex{4})
        ec = default_embedded_cubature(reference_domain(domain))
        @test typeof(ec) <: EmbeddedCubature
    end
end

@testset "Integrate" begin
    domain = reference_domain(Segment)

    buffer = allocate_buffer(x -> sum(x), domain)
    @test typeof(buffer) <: BinaryHeap

    I, E = integrate(x -> sin(10 * x[1]), domain; buffer=buffer)
    R = sin(5)^2 / 5
    @test abs(I - R) ≤ E * abs(R)

    I, E = integrate(x -> cos(7.5 * x[1]), domain; buffer=buffer)
    R = 2 * sin(7.5) / 15
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a segment" begin
    domain = Segment{float(Int)}(0, 1)
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (
        embedded_cubature(SEGMENT_GK7),
        default_embedded_cubature(domain),
        embedded_cubature(SEGMENT_GK31),
    )
        for (fct, R) in [
            (x -> exp(x[1]), exp(1) - 1),
            (x -> cos(10 * x[1]), sin(10) / 10),
            (x -> 1 / √x[1], 2),
        ]
            I, E = integrate(fct, domain; buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a triangle" begin
    domain = triangle((0, 0), (2, 0), (0, 2))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (embedded_cubature(TRIANGLE_GM19), default_embedded_cubature(domain))
        for (fct, R) in [
            (x -> exp(x[1] + 3 * x[2]), (exp(6) - 3 * exp(2) + 2) / 6),
            (x -> cos(7 * x[1] + 3 * x[2]), (-3 * cos(14) + 7 * cos(6) - 4) / 84),
            (x -> 1 / norm(x), 2 * sqrt(2) * asinh(1)),
        ]
            I, E = integrate(fct, domain; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a rectangle" begin
    domain = rectangle((0, 0), (1, 1))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (
        embedded_cubature(SQUARE_GM17),
        embedded_cubature(SQUARE_CH21),
        default_embedded_cubature(domain),
    )
        for (fct, R) in [
            (x -> exp(x[1] + x[2]), (exp(1) - 1)^2),
            (x -> 1 / norm(x), log(17 + 12 * sqrt(2)) / 2),
        ]
            I, E = integrate(fct, domain; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a tetrahedron" begin
    domain = tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for (fct, R) in [
        (
            x -> exp(x[1] + 3 * x[2] + 5 * x[3]),
            (3 * exp(5) - 10 * exp(3) + 15 * exp(1) - 8) / 120,
        ),
        (x -> 1 / norm(x), 0.3614258523411),
    ]
        I, E = integrate(fct, domain; buffer=buffer)
        @test abs(I - R) ≤ E * abs(R)
    end
end

@testset "Integrate over a Cuboid" begin
    domain = cuboid((0, 0, 0), (1, 1, 1))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (
        embedded_cubature(CUBE_GM33),
        default_embedded_cubature(domain),
        embedded_cubature(CUBE_BE115),
    )
        for (fct, R) in [
            (x -> exp(x[1]), exp(1) - 1),
            (x -> 1 / (1 + norm(x)^2)^2, π^2 / 32),
            (x -> 1 / norm(x), 1.1900386819897766),
        ]
            I, E = integrate(fct, domain; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "4-Simplex" begin
    I, E = integrate(x -> 1 / norm(x), reference_domain(Simplex{4}))
    R = 0.089876019011
    @test abs(I - R) ≤ E * abs(R)
end

@testset "4-Orthotope" begin
    I, E = integrate(x -> 1 / norm(x), reference_domain(Orthotope{4}))
    R = 0.9674120212411487
    @test abs(I - R) ≤ E * abs(R)
end
