using ForwardDiff
using HAdaptiveIntegration
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using LinearAlgebra
using Logging
using Test

global_logger(SimpleLogger(stderr, Logging.Warn))

@testset "IncreasePrecisionExt" begin
    @testset "Triangle" begin
        ec_ref = embedded_cubature(TRIANGLE_RL19, Float64)

        ec_f32 = embedded_cubature(TRIANGLE_RL19, Float32)
        tec0 = TabulatedEmbeddedCubature{Triangle}(;
            description="Reduce Radon-Laurie",
            reference="",
            order_high=TRIANGLE_RL19.order_high,
            order_low=TRIANGLE_RL19.order_low,
            precision=7,
            nodes=[Vector([string(v) for v in p]) for p in ec_f32.nodes],
            weights_high=string.(ec_f32.weights_high),
            weights_low=string.(ec_f32.weights_low),
        )
        ec0 = embedded_cubature(tec0, Float64)

        tec1 = increase_precision(tec0, Float64; atol=1e-14)
        @test typeof(tec1) <: TabulatedEmbeddedCubature{Triangle}
        ec1 = embedded_cubature(tec1, Float64)

        domain = reference_domain(Triangle{Float64})
        fct = x -> exp(2 * x[1] + x[2])
        R = (exp(1) - 1)^2 / 2
        rtol = 1e-12

        I, E = integrate(fct, domain; embedded_cubature=ec0, rtol=rtol)
        @show I, E, R
        @test abs(I - R) > rtol * abs(R)

        I, E = integrate(fct, domain; embedded_cubature=ec1, rtol=rtol)
        @show I, E, R
        @test abs(I - R) ≤ rtol * abs(R)
    end

    @testset "Rectangle" begin
        tec0 = simple_cbt_orthotope(2, Float32)

        tec1 = increase_precision(tec0, Float64)
        @test typeof(tec1) <: TabulatedEmbeddedCubature{Orthotope{2}}

        ec = embedded_cubature(tec1, Float64)
    end
end
