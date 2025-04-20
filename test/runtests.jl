using Aqua
using HAdaptiveIntegration
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HAdaptiveIntegration)
end

@testset "HAdaptiveIntegration.jl" begin
    include("test_domain.jl")
    include("test_rule.jl")
    include("test_integrate.jl")
end
