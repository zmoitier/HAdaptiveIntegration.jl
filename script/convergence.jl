using CairoMakie
using HAdaptiveIntegration
using LinearAlgebra: norm

function main()
    triangle = Triangle((0, 0), (1, 0), (0, 1))
    fct = x -> 1 / norm(x)
    Iₑₓ = sqrt(2) * asinh(1)

    rule = HAdaptiveIntegration.Rule.RadonLaurie()

    ec = HAdaptiveIntegration.embedded_cubature(rule)
    nb_nodes = length(ec.nodes)

    nb_evaluations = Vector{Float64}()
    estimated_values = Vector{Float64}()
    estimated_errors = Vector{Float64}()

    function callback(I, E, nb_subdiv, buffer)
        push!(nb_evaluations, (nb_subdiv + 1) * nb_nodes)
        push!(estimated_values, I)
        push!(estimated_errors, E)
        return nothing
    end

    I, E = integrate(fct, triangle; embedded_cubature=ec, callback=callback)
    absolute_error = abs.(estimated_values .- Iₑₓ)

    k = round(Int, 0.25 * length(nb_evaluations))
    @show slope_I = log10.(nb_evaluations[k:end]) \ log10.(absolute_error[k:end])
    @show slope_E = log10.(nb_evaluations[k:end]) \ log10.(estimated_errors[k:end])

    fig = Figure()
    ax = Axis(
        fig[1, 1]; xlabel="Number of function evaluations", xscale=log10, yscale=log10
    )

    scatterlines!(ax, nb_evaluations, absolute_error; label="|I - Iₑₓ|")
    scatterlines!(ax, nb_evaluations, estimated_errors; label="E")

    b = last(nb_evaluations)
    lines!(
        ax,
        nb_evaluations[k:end],
        last(estimated_errors) * (nb_evaluations[k:end] ./ b) .^ slope_E;
        color=:black,
        label="order $slope_E",
    )
    lines!(
        ax,
        nb_evaluations[k:end],
        last(absolute_error) * (nb_evaluations[k:end] ./ b) .^ slope_I;
        linestyle=:dash,
        color=:black,
        label="order $slope_I",
    )

    axislegend(ax)

    return fig
end

main()
