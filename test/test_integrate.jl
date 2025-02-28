using Test, StaticArrays
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
    # Test on Float64 precision
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
    # Test on Float64 precision
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
