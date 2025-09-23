using DataStructures: BinaryHeap
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using HAdaptiveIntegration:
    allocate_buffer, default_embedded_cubature, default_subdivision, integrate, resum
using LinearAlgebra: norm
using Test

@testset "Default" begin
    for domain_type in
        (Segment, Simplex{1}, Triangle, Tetrahedron, Orthotope{1}, Rectangle, Cuboid)
        domain = reference_domain(domain_type)

        @test typeof(default_subdivision(domain)) <: Function
        @test typeof(default_embedded_cubature(domain)) <: EmbeddedCubature
    end
end

@testset "Integrate over a segment" begin
    domain = Segment(0, 1)
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (
        embedded_cubature(SEGMENT_GK7),
        default_embedded_cubature(domain),
        embedded_cubature(SEGMENT_GK31),
    )
        for (fct, R) in [(x -> exp(x[1]), exp(1) - 1), (x -> 1 / sqrt(x[1]), 2)]
            I, E = integrate(fct, domain; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a triangle" begin
    domain = Triangle((0, 0), (2, 0), (0, 2))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (embedded_cubature(TRIANGLE_GM19), default_embedded_cubature(domain))
        for (fct, R) in [
            (x -> cos(7 * x[1] + 3 * x[2]), (-3 * cos(14) + 7 * cos(6) - 4) / 84),
            (x -> 1 / norm(x), 2 * sqrt(2) * asinh(1)),
        ]
            I, E = integrate(fct, domain; embedded_cubature=ec, buffer=buffer)
            @test abs(I - R) ≤ E * abs(R)
        end
    end
end

@testset "Integrate over a rectangle" begin
    domain = Rectangle((0, 0), (1, 1))
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
    domain = Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
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
    domain = Cuboid((0, 0, 0), (1, 1, 1))
    buffer = allocate_buffer(x -> zero(x[1]), domain)

    for ec in (
        embedded_cubature(CUBE_GM33),
        default_embedded_cubature(domain),
        embedded_cubature(CUBE_BE115),
    )
        for (fct, R) in
            [(x -> 1 / (1 + norm(x)^2)^2, π^2 / 32), (x -> 1 / norm(x), 1.1900386819897766)]
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

@testset "Return buffer" begin
    domain = Segment(0, 1)
    fct = x -> exp(x[1])

    # by default, the buffer is not returned
    @test length(integrate(fct, domain)) == 2

    # we can ask for the buffer to be returned
    @test length(integrate(fct, domain; return_buffer=Val(true))) == 3
    I, E, buf = integrate(fct, domain; return_buffer=Val(true))
    @test typeof(buf) == typeof(allocate_buffer(fct, domain))

    I2, E2 = resum(buf)
    @test I ≈ I2
    @test E ≈ E2
end
