using BenchmarkTools
using HAdaptiveIntegration
using HCubature
using LinearAlgebra
using StaticArrays

function uniform_regularity(dim::Int)
    fct = x -> 1 + prod(cos.(4 .* x .+ 1))

    a, b = zeros(dim), ones(dim)

    fct_count = x -> (counter[] += 1; fct(x))

    counter = Ref(0)
    I, E = integrate(fct_count, Orthotope(a, b); rtol=1e-8)
    @show I E counter[]

    counter[] = 0
    I, E = hcubature(fct_count, a, b; rtol=1e-8)
    @show I E counter[]

    return nothing
end

function singularity_point(dim::Int)
    x₀ = SVector{dim}([-1e-2, zeros(dim - 1)...])

    fct = x -> cis(4 * norm(x - x₀)) / norm(x - x₀)
    a, b = zeros(dim), ones(dim)

    fct_count = x -> (counter[] += 1; fct(x))

    counter = Ref(0)
    I, E = integrate(fct_count, Orthotope(a, b); rtol=1e-8)
    @show I E counter[]

    counter[] = 0
    I, E = hcubature(fct_count, a, b; rtol=1e-8)
    @show I E counter[]

    return nothing
end

function singularity_hyperplane(dim::Int)
    x₀ = SVector{dim - 1}([-1e-2, zeros(dim - 2)...])

    fct = x -> cis(4 * norm(x[1:dim-1] - x₀)) / norm(x[1:dim-1] - x₀)
    a, b = zeros(dim), ones(dim)

    fct_count = x -> (counter[] += 1; fct(x))

    counter = Ref(0)
    I, E = integrate(fct_count, Orthotope(a, b); rtol=1e-8)
    @show I E counter[]

    counter[] = 0
    I, E = hcubature(fct_count, a, b; rtol=1e-8)
    @show I E counter[]

    return nothing
end
