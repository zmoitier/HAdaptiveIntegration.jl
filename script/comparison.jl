using BenchmarkTools
using HAdaptiveIntegration: HAdaptiveIntegration as HAI
using HCubature: HCubature as HC
using LinearAlgebra

function square()
    fct = x -> 1 + prod(cos.(32 .* x))
    # fct = x -> sqrt(norm(x))
    # fct = x -> sqrt(norm(x[1]))

    a = (0.0, 0.0)
    b = (1.0, 1.0)
    dom = HAI.Rectangle(a, b)

    εᵣₑₗ = sqrt(eps(Float64))

    @show HAI.integrate(fct, dom; rtol=εᵣₑₗ)
    display(@benchmark HAI.integrate($fct, $dom; rtol=$εᵣₑₗ))

    println()

    @show HC.hcubature(fct, a, b; rtol=εᵣₑₗ)
    display(@benchmark HC.hcubature($fct, $a, $b; rtol=$εᵣₑₗ))

    return nothing
end

function cube()
    fct = x -> 1 + prod(cos.(32 .* x))
    # fct = x -> sqrt(norm(x))
    # fct = x -> sqrt(norm(x[1]))
    # fct = x -> sqrt(norm(x[1:2]))

    a = (0.0, 0.0, 0.0)
    b = (1.0, 1.0, 1.0)
    dom = HAI.Cuboid(a, b)

    εᵣₑₗ = sqrt(eps(Float64))

    @show HAI.integrate(fct, dom; rtol=εᵣₑₗ)
    display(@benchmark HAI.integrate($fct, $dom; rtol=$εᵣₑₗ))

    println()

    @show HC.hcubature(fct, a, b; rtol=εᵣₑₗ)
    display(@benchmark HC.hcubature($fct, $a, $b; rtol=$εᵣₑₗ))

    return nothing
end

function tesseract()
    fct = x -> 1 + prod(cos.(32 .* x))
    # fct = x -> sqrt(norm(x))
    # fct = x -> sqrt(norm(x[1]))
    # fct = x -> sqrt(norm(x[1:2]))
    # fct = x -> sqrt(norm(x[1:3]))

    a = (0.0, 0.0, 0.0, 0.0)
    b = (1.0, 1.0, 1.0, 1.0)
    dom = HAI.Orthotope(a, b)

    εᵣₑₗ = eps(Float64)^0.25

    @show HAI.integrate(fct, dom; rtol=εᵣₑₗ)
    display(@benchmark HAI.integrate($fct, $dom; rtol=$εᵣₑₗ))

    println()

    @show HC.hcubature(fct, a, b; rtol=εᵣₑₗ)
    display(@benchmark HC.hcubature($fct, $a, $b; rtol=$εᵣₑₗ))

    return nothing
end
