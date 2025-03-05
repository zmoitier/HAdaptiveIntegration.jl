using Test, StaticArrays, LinearAlgebra
import HAdaptiveIntegration as HAI

@testset "Integrate over a segment" begin
    # Test on Float64 precision
    @test isnothing(HAI.check_order(HAI.SEGMENT_G7K15, HAI.Segment{Float64}))

    ec = HAI.embedded_cubature(HAI.SEGMENT_G7K15, Float64)
    segment = HAI.segment(0.0, 1.0)

    I, E = ec(x -> exp(x[1]), segment)
    R = exp(1) - exp(0)
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> cos(10 * x[1]), segment)
    R = sin(10) / 10
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> 1 / √x[1], segment)
    R = 2
    @test abs(I - R) ≤ E * abs(R)

    I, E = HAI.integrate(x -> 1 / √x[1], segment, ec)
    R = 2
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a rectangle" begin
    @test isnothing(HAI.check_order(HAI.SQUARE_CH21G25, HAI.Rectangle{Float64}))

    ec = HAI.embedded_cubature(HAI.SQUARE_CH21G25, Float64)
    square = HAI.rectangle((0.0, 0.0), (1.0, 1.0))

    # FIXME: it seems the Cools Haegemens rule underestimates the error for one-dimensional integrands
    I, E = ec(x -> exp(x[1]), square)
    R = exp(1) - exp(0)
    @test abs(I - R) < 1e-8
    @test_broken abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> exp(x[2]), square)
    R = exp(1) - exp(0)
    @test abs(I - R) < 1e-8
    @test_broken abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> cos(10 * x[1]), square)
    R = sin(10) / 10
    @test abs(I - R) ≤ 1e-3
    @test_broken abs(I - R) ≤ E * abs(R)

    # Works for two-dimensional integrands
    I, E = ec(x -> exp(x[1]) * exp(x[2]), square)
    R = (exp(1) - exp(0))^2
    @test abs(I - R) ≤ E * abs(R)

    I, E = ec(x -> 1 / norm(x), square)
    R = 1 / 2 * log(17 + 12 * sqrt(2))
    @test abs(I - R) ≤ E * abs(R)

    I, E = HAI.integrate(x -> 1 / norm(x), square, ec)
    @test abs(I - R) < 1e-8
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a triangle" begin
    @test isnothing(HAI.check_order(HAI.TRIANGLE_R7L19, HAI.Triangle{Float64}))

    ec = HAI.embedded_cubature(HAI.TRIANGLE_R7L19, Float64)

    triangle = HAI.triangle((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
    I, E = ec(x -> exp(x[1] + 3 * x[2]), triangle)
    R = (exp(6) - 3 * exp(2) + 2) / 6
    @test abs(I - R) ≤ E * abs(R)

    triangle = HAI.triangle((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
    I, E = ec(x -> cos(7 * x[1] + 3 * x[2]), triangle)
    R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
    @test abs(I - R) ≤ E * abs(R)

    triangle = HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    I, E = HAI.integrate(x -> 1 / norm(x), triangle, ec)
    R = sqrt(2) * asinh(1)
    @test abs(I - R) ≤ E * abs(R)
end

@testset "Integrate over a Cuboid" begin
    @test isnothing(HAI.check_order(HAI.CUBE_BE65, HAI.Cuboid{Float64}))

    ec = HAI.embedded_cubature(HAI.CUBE_BE65, Float64)
    cube = HAI.cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))

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

    I, E = HAI.integrate(x -> 1 / norm(x), cube, ec)
    @test E ≤ √eps(Float64) * I
end

@testset "GrundmannMoeller quadrature" begin
    @testset "Triangle" begin
        degree = 7 # must be odd
        ec = HAI.embedded_cubature(HAI.GrundmannMoeller(), 2, degree, Float64)
        @test isnothing(
            HAI.check_order(
                ec, degree, degree - 2, HAI.Triangle{Float64}; rtol=15 * eps(Float64)
            ),
        )

        rtol = 1e-8

        triangle = HAI.triangle((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
        I, E = HAI.integrate(x -> exp(x[1] + 3 * x[2]), triangle, ec; rtol)
        R = (exp(6) - 3 * exp(2) + 2) / 6
        @test abs(I - R) ≤ rtol * abs(R)

        triangle = HAI.triangle((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
        I, E = HAI.integrate(x -> cos(7 * x[1] + 3 * x[2]), triangle, ec; rtol)
        R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
        @test abs(I - R) ≤ rtol * abs(R)

        triangle = HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        I, E = HAI.integrate(x -> 1 / norm(x), triangle, ec; rtol)
        R = sqrt(2) * asinh(1)
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "Tetrahedron" begin
        degree = 7 # must be odd
        ec = HAI.embedded_cubature(HAI.GrundmannMoeller(), 3, degree, Float64)
        @test isnothing(
            HAI.check_order(
                ec, degree, degree - 2, HAI.Tetrahedron{Float64}; rtol=15 * eps(Float64)
            ),
        )

        rtol = 1e-8
        tetrahedron = HAI.tetrahedron(
            (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)
        )

        I, E = HAI.integrate(x -> 1, tetrahedron, ec; rtol)
        R = 1 / 6
        @test abs(I - R) ≤ rtol * abs(R)

        I, E = HAI.integrate(x -> exp(x[1] + 3 * x[2] + 5 * x[3]), tetrahedron, ec; rtol)
        R = (3 * exp(5) - 10 * exp(3) + 15 * exp(1) - 8) / 120
        @test abs(I - R) ≤ rtol * abs(R)

        I, E = HAI.integrate(x -> 1 / norm(x), tetrahedron, ec; rtol)
        R = 0.361426 # not sure if an analytic solution exists... this was from WolframAlpha
        @test abs(I - R) ≤ 1e-4
    end

    @testset "4-Simplex" begin
        degree = 7 # must be odd
        ec = HAI.embedded_cubature(HAI.GrundmannMoeller(), 4, degree, Float64)
        @test isnothing(
            HAI.check_order(
                ec, degree, degree - 2, HAI.Simplex{4,Float64,5}; rtol=50 * eps(Float64)
            ),
        )
    end
end
