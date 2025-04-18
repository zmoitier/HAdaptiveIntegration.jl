using Aqua
using HAdaptiveIntegration
using Test

@testset "Code quality (Aqua.jl)" begin
    Aqua.test_all(HAdaptiveIntegration)
end
