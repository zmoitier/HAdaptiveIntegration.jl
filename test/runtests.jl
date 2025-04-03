using Test

@testset "HAdaptiveIntegration.jl" begin
    include("aqua_test.jl")

    include("test_domain.jl")
    include("test_rule.jl")
    include("test_integrate.jl")
end
