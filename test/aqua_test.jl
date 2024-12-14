using HAdaptiveIntegration
using Test
using Aqua

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HAdaptiveIntegration)
end
