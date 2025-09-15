using ForwardDiff
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule
using HAdaptiveIntegration: increase_precision
using Logging
using Test

global_logger(SimpleLogger(stderr, Logging.Warn))

@testset "IncreasePrecisionExt" begin
    @testset "Simplex" begin
        tec = TabulatedEmbeddedCubature{Triangle}(;
            description="Custom with 4 nodes",
            reference="direct computation",
            order_high=2,
            order_low=1,
            precision=6,
            nodes=[
                ["$(Float32(1/3))", "$(Float32(1/3))"],
                ["$(Float32(0))", "$(Float32(0))"],
                ["$(Float32(1))", "$(Float32(0))"],
                ["$(Float32(0))", "$(Float32(1))"],
            ],
            weights_high=[
                "$(Float32(3/8))",
                "$(Float32(1/24))",
                "$(Float32(1/24))",
                "$(Float32(1/24))",
            ],
            weights_low=["$(Float32(1/2))"],
        )

        tec_incr = increase_precision(tec, Float64)
        @test typeof(tec_incr) <: TabulatedEmbeddedCubature{Triangle}

        ec = embedded_cubature(tec_incr, Float64)
    end

    @testset "Orthotope" begin
        tec = TabulatedEmbeddedCubature{Rectangle}(;
            description="Custom with 5 nodes",
            reference="direct computation",
            order_high=3,
            order_low=1,
            precision=6,
            nodes=[
                ["$(Float32(1/2))", "$(Float32(1/2))"],
                ["$(Float32(0))", "$(Float32(0))"],
                ["$(Float32(1))", "$(Float32(0))"],
                ["$(Float32(0))", "$(Float32(1))"],
                ["$(Float32(1))", "$(Float32(1))"],
            ],
            weights_high=[
                "$(Float32(2/3))",
                "$(Float32(1/12))",
                "$(Float32(1/12))",
                "$(Float32(1/12))",
                "$(Float32(1/12))",
            ],
            weights_low=["$(Float32(1))"],
        )

        tec_incr = increase_precision(tec, Float64)
        @test typeof(tec_incr) <: TabulatedEmbeddedCubature{Rectangle}

        ec = embedded_cubature(tec_incr, Float64)
    end
end
