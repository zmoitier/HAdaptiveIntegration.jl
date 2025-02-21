using Test, LinearAlgebra

@testset "HAdaptiveIntegration.jl" begin
    # include("aqua_test.jl")

    include("test_simplex.jl")
    include("test_orthotope.jl")

    @testset "Segment" begin
        segment = HAI.Segment((0.0,), (1.0,))
        # Test on Float64 precision
        T = Float64
        Q = HAI.Quadrature{T}
        @test HAI.check_order(Q(HAI.SEGMENT_GAUSS_O13_N7), 13)
        @test HAI.check_order(Q(HAI.SEGMENT_KRONROD_O23_N15), 23)
        embd_quad = HAI.EmbeddedQuadrature(; name="segment-G7K15", datatype=T)

        I, E = embd_quad(x -> exp(x[1]), segment)
        R = exp(1) - exp(0)
        @test abs(I - R) ≤ E * abs(R)

        I, E = embd_quad(x -> cos(10 * x[1]), segment)
        R = sin(10) / 10
        @test abs(I - R) ≤ E * abs(R)

        I, E = embd_quad(x -> 1 / √x[1], segment)
        R = 2
        @test abs(I - R) ≤ E * abs(R)

        I, E = HAI.integrate(x -> 1 / √x[1], segment, embd_quad)
        R = 2
        @test abs(I - R) ≤ E * abs(R)
    end

    @testset "Triangle" begin
        triangle = HAI.Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

        # Test on Float64 precision
        T = Float64
        Q = HAI.Quadrature{T}
        @test HAI.check_order(Q(HAI.TRIANGLE_RADON_O5_N7), 5)
        @test HAI.check_order(Q(HAI.TRIANGLE_LAURIE_O8_N19), 8)

        embd_quad = HAI.EmbeddedQuadrature(; name="triangle-LaurieRadon", datatype=T)

        simplex = HAI.Simplex((0, 0), (2, 0), (0, 2))
        I, E = embd_quad(x -> exp(x[1] + 3 * x[2]), simplex)
        R = (exp(6) - 3 * exp(2) + 2) / 6
        @test abs(I - R) ≤ E * abs(R)

        simplex = HAI.Simplex((0, 0), (2, 0), (0, 2))
        I, E = embd_quad(x -> cos(7 * x[1] + 3 * x[2]), simplex)
        R = (-3 * cos(14) + 7 * cos(6) - 4) / 84
        @test abs(I - R) ≤ E * abs(R)

        simplex = HAI.Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
        I, E = HAI.integrate(x -> 1 / norm(x), simplex, embd_quad)
        R = sqrt(2) * asinh(1)
        @test abs(I - R) ≤ E * abs(R)
    end

    @testset "Square" begin
        sq = HAI.Rectangle((0.0, 0.0), (1.0, 1.0))

        # Test on Float64 precision
        T = Float64
        Q1 = HAI.Quadrature{T}(HAI.SQUARE_GAUSS_O9_N25)
        @test HAI.check_order_square(Q1, 9)
        Q2 = HAI.Quadrature{T}(HAI.SQUARE_COOLS_HAEGEMANS_O7_N21)
        @test HAI.check_order_square(Q2, 7)

        embd_quad = HAI.EmbeddedQuadrature(; name="square-CoolsHaegemans", datatype=T)

        # FIXME: it seems the CoolsHaegemens rule underestimates the
        # error for one-dimensional integrands
        I, E = embd_quad(x -> exp(x[1]), sq)
        R = exp(1) - exp(0)
        @test abs(I - R) < 1e-8
        @test_broken abs(I - R) ≤ E * abs(R)
        I, E = embd_quad(x -> exp(x[2]), sq)
        R = exp(1) - exp(0)
        @test abs(I - R) < 1e-8
        @test_broken abs(I - R) ≤ E * abs(R)

        I, E = embd_quad(x -> cos(10 * x[1]), sq)
        R = sin(10) / 10
        @test abs(I - R) ≤ 1e-3
        @test_broken abs(I - R) ≤ E * abs(R)

        # Works for two-dimensional integrands
        I, E = embd_quad(x -> exp(x[1]) * exp(x[2]), sq)
        R = (exp(1) - exp(0))^2
        @test abs(I - R) ≤ E * abs(R)

        I, E = embd_quad(x -> 1 / norm(x), sq)
        R = 1 / 2 * log(17 + 12 * sqrt(2))
        @test abs(I - R) ≤ E * abs(R)
        I, E = HAI.integrate(x -> 1 / norm(x), sq, embd_quad)
        @test abs(I - R) < 1e-8
        @test abs(I - R) ≤ E * abs(R)
    end
end
