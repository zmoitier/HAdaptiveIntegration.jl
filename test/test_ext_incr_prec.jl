using ForwardDiff
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using HAdaptiveIntegration: increase_precision, integrate
using LinearAlgebra
using Logging
using Test

global_logger(SimpleLogger(stderr, Logging.Warn))

@testset "IncreasePrecisionExt" begin
    @testset "Segment" begin
        domain = reference_domain(Segment{Float64})
        fct = x -> exp(x[1])
        R = exp(1) - 1
        rtol = 1e-14

        ec_f32 = embedded_cubature(SEGMENT_GK15, Float32)
        tec0 = TabulatedEmbeddedCubature{Segment}(;
            description="Reduce Gauss-Kronrod",
            reference="",
            order_high=SEGMENT_GK15.order_high,
            order_low=SEGMENT_GK15.order_low,
            precision=7,
            nodes=[Vector([string(v) for v in p]) for p in ec_f32.nodes],
            weights_high=string.(ec_f32.weights_high),
            weights_low=string.(ec_f32.weights_low),
        )

        ec0 = embedded_cubature(tec0, Float64)
        I, E = integrate(fct, domain; embedded_cubature=ec0, rtol=rtol)
        @test abs(I - R) > rtol * abs(R)

        setprecision(BigFloat, 32; base=10)
        tec1 = increase_precision(tec0, BigFloat)
        @test typeof(tec1) <: TabulatedEmbeddedCubature{Segment}
        ec1 = embedded_cubature(tec1, Float64)

        I, E = integrate(fct, domain; embedded_cubature=ec1, rtol=rtol)
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "Triangle" begin
        domain = reference_domain(Triangle{Float64})
        fct = x -> exp(2 * x[1] + x[2])
        R = (exp(1) - 1)^2 / 2
        rtol = 1e-14

        ec_f32 = embedded_cubature(RadonLaurie(), Float32)
        oh, ol = orders(RadonLaurie())
        tec0 = TabulatedEmbeddedCubature{Triangle}(;
            description="Reduce Radon-Laurie",
            reference="",
            order_high=oh,
            order_low=ol,
            precision=7,
            nodes=[Vector([string(v) for v in p]) for p in ec_f32.nodes],
            weights_high=string.(ec_f32.weights_high),
            weights_low=string.(ec_f32.weights_low),
        )

        ec0 = embedded_cubature(tec0, Float64)
        I, E = integrate(fct, domain; embedded_cubature=ec0, rtol=rtol)
        @test abs(I - R) > rtol * abs(R)

        setprecision(BigFloat, 20; base=10)
        tec1 = increase_precision(tec0, BigFloat)
        @test typeof(tec1) <: TabulatedEmbeddedCubature{Triangle}
        ec1 = embedded_cubature(tec1, Float64)

        I, E = integrate(fct, domain; embedded_cubature=ec1, rtol=rtol)
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "Rectangle" begin
        domain = reference_domain(Rectangle{Float64})
        fct = x -> exp(2 * x[1] + x[2])
        R = (exp(2) - 1) * (exp(1) - 1) / 2
        rtol = 1e-14

        ec_f32 = embedded_cubature(SQUARE_CH25, Float32)
        tec0 = TabulatedEmbeddedCubature{Rectangle}(;
            description="Reduce Cools-Haegemans",
            reference="",
            order_high=SQUARE_CH25.order_high,
            order_low=SQUARE_CH25.order_low,
            precision=7,
            nodes=[Vector([string(v) for v in p]) for p in ec_f32.nodes],
            weights_high=string.(ec_f32.weights_high),
            weights_low=string.(ec_f32.weights_low),
        )

        ec0 = embedded_cubature(tec0, Float64)
        I, E = integrate(fct, domain; embedded_cubature=ec0, rtol=rtol)
        @test abs(I - R) > rtol * abs(R)

        setprecision(BigFloat, 64; base=10)
        tec1 = increase_precision(tec0, BigFloat; f_atol=1e-25)
        @test typeof(tec1) <: TabulatedEmbeddedCubature{Rectangle}
        ec1 = embedded_cubature(tec1, Float64)

        I, E = integrate(fct, domain; embedded_cubature=ec1, rtol=rtol)
        @test abs(I - R) ≤ rtol * abs(R)
    end
end
>
