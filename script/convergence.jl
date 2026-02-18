using CairoMakie
using HAdaptiveIntegration
using LinearAlgebra
using Printf

function log_linear_regression(x::Vector{T}, y::Vector{T}) where {T<:Real}
    A = hcat(ones(T, length(x)), log10.(x))
    c = A \ log10.(y)
    return (10^c[1], c[2])
end

function main()
    triangle = Triangle((0, 0), (2, 0), (0, 2))

    fct = x -> 1 / norm(x)
    Iₑₓ = 2 * sqrt(2) * asinh(1)

    # fct = x -> cos(7 * x[1] + 3 * x[2])
    # Iₑₓ = (-3 * cos(14) + 7 * cos(6) - 4) / 84

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
    absolute_error = abs.(estimated_values .- Iₑₓ) .+ 1e-16

    k = round(Int, 0.6 * length(nb_evaluations))
    c_I, s_I = log_linear_regression(nb_evaluations[k:end], absolute_error[k:end])
    c_E, s_E = log_linear_regression(nb_evaluations[k:end], estimated_errors[k:end])

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
        (√10 * c_I) * nb_evaluations[k:end] .^ s_I;
        color=:black,
        label=@sprintf("slope %.2f", s_I),
    )
    lines!(
        ax,
        nb_evaluations[k:end],
        (√10 * c_E) * nb_evaluations[k:end] .^ s_E;
        linestyle=:dash,
        color=:black,
        label=@sprintf("slope %.2f", s_E),
    )

    axislegend(ax)

    return fig
end

main()
