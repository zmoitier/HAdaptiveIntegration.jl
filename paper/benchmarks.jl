using BenchmarkTools
using CairoMakie
using HAdaptiveIntegration
using HAdaptiveIntegration.Rule: CUBE_BE65, GenzMalik, SQUARE_CH25, embedded_cubature
using HCubature
using LaTeXStrings
using LinearAlgebra
using Quadmath
using StaticArrays

function peak(x, c)
    return c * exp(-c * x^2)
end

function reference(f::FCT, ::Val{D}) where {FCT, D}
    rule = if D == 2
        embedded_cubature(SQUARE_CH25, Float128)
    elseif D == 3
        embedded_cubature(CUBE_BE65, Float128)
    else
        embedded_cubature(GenzMalik{D}(), Float128)
    end
    a, b = zeros(Float128, D), ones(Float128, D)
    I, E = integrate(f, Orthotope(a, b); rtol = Float128(1.0e-14), rule = rule)
    return Float64(I), Float64(E)
end

function plot_nb_eval(rtol_vec, hai, hc)
    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xscale = log10,
        xlabel = "rtol",
        yscale = log10,
        ylabel = "Evaluation count",
    )

    scatterlines!(ax, rtol_vec, hai.N; marker = :cross, linestyle = :dash, label = "HAI")
    scatterlines!(ax, rtol_vec, hc.N; marker = :xcross, linestyle = :dash, label = "HC")

    axislegend(ax; position = :rt)

    return fig
end

function plot_error(rtol_vec, hai, hc, Iᵣ)
    colors = Makie.wong_colors()

    fig = Figure()
    ax = Axis(
        fig[1, 1];
        xscale = log10,
        xlabel = "rtol",
        yscale = log10,
        ylabel = "Relative error",
    )

    scatterlines!(
        ax,
        rtol_vec,
        abs.(hai.I .- Iᵣ) ./ abs(Iᵣ);
        marker = :cross,
        linestyle = :dash,
        color = colors[1],
        label = "HAI true error",
    )
    scatterlines!(
        ax,
        rtol_vec,
        abs.(hc.I .- Iᵣ) ./ abs(Iᵣ);
        marker = :xcross,
        linestyle = :dash,
        color = colors[2],
        label = "HC true error"
    )

    scatterlines!(
        ax,
        rtol_vec,
        hai.E ./ abs(Iᵣ),
        marker = :cross,
        linestyle = :dot,
        color = colors[1],
        label = "HAI: estimated error",
    )
    scatterlines!(
        ax,
        rtol_vec,
        hc.E ./ abs(Iᵣ),
        marker = :xcross,
        linestyle = :dot,
        color = colors[2],
        label = "HC: estimated error",
    )

    axislegend(ax; position = :lt)

    return fig
end

function uniform_subdivision(dim::Int)
    fct = if dim == 2
        x -> 1 + prod(cos.(8 .* x .+ 1))
    elseif dim == 3
        x -> 1 + prod(cos.(4 .* x .+ 1))
    else
        x -> 1 + prod(cos.(2 .* x .+ 1))
    end

    a, b = zeros(dim), ones(dim)

    Iᵣ, Eᵣ = reference(fct, Val(dim))
    println("Reference: I = $Iᵣ, E = $Eᵣ")

    rtol_vec = 10.0 .^ range(-1, stop = -10, step = -1)

    hai = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )
    hc = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )

    counter = Ref(0)
    fct_count = x -> (counter[] += 1; fct(x))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, Orthotope(a, b); rtol = rtol)
        hai.I[i] = I
        hai.E[i] = E
        hai.N[i] = counter[]

        counter[] = 0
        I, E = hcubature(fct_count, a, b; rtol = rtol)
        hc.I[i] = I
        hc.E[i] = E
        hc.N[i] = counter[]
    end

    plot_nb_eval(rtol_vec, hai, hc) |> display
    plot_error(rtol_vec, hai, hc, Iᵣ) |> display

    return nothing
end

function singularity_point(dim::Int)
    x₀ = SVector{dim}(1 ./ range(3, length = dim, step = 2))
    fct = if dim == 2
        x -> peak(norm(x - x₀), 8)
    elseif dim == 3
        x -> peak(norm(x - x₀), 4)
    else
        x -> peak(norm(x - x₀), 2)
    end

    a, b = zeros(dim), ones(dim)

    Iᵣ, Eᵣ = reference(fct, Val(dim))
    println("Reference: I = $Iᵣ, E = $Eᵣ")

    rtol_vec = 10.0 .^ range(-1, stop = -10, step = -1)

    hai = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )
    hc = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )

    counter = Ref(0)
    fct_count = x -> (counter[] += 1; fct(x))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, Orthotope(a, b); rtol = rtol)
        hai.I[i] = I
        hai.E[i] = E
        hai.N[i] = counter[]

        counter[] = 0
        I, E = hcubature(fct_count, a, b; rtol = rtol)
        hc.I[i] = I
        hc.E[i] = E
        hc.N[i] = counter[]
    end

    plot_nb_eval(rtol_vec, hai, hc) |> display
    plot_error(rtol_vec, hai, hc, Iᵣ) |> display

    return nothing
end

function singularity_hyperplane(dim::Int)
    x₀ = 1 / 3 # SVector{dim}(1 ./ range(3, length = dim, step = 2))

    fct = if dim == 2
        x -> (1 + prod(cos.(2 .* x .+ 1))) * peak(x[1] - x₀, 16)
    elseif dim == 3
        x -> (1 + prod(cos.(2 .* x .+ 1))) * peak(x[1] - x₀, 8)
    else
        x -> (1 + prod(cos.(2 .* x .+ 1))) * peak(x[1] - x₀, 4)
    end

    a, b = zeros(dim), ones(dim)

    Iᵣ, Eᵣ = reference(fct, Val(dim))
    println("Reference: I = $Iᵣ, E = $Eᵣ")

    rtol_vec = 10.0 .^ range(-1, stop = -10, step = -1)

    hai = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )
    hc = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )

    counter = Ref(0)
    fct_count = x -> (counter[] += 1; fct(x))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, Orthotope(a, b); rtol = rtol)
        hai.I[i] = I
        hai.E[i] = E
        hai.N[i] = counter[]

        counter[] = 0
        I, E = hcubature(fct_count, a, b; rtol = rtol)
        hc.I[i] = I
        hc.E[i] = E
        hc.N[i] = counter[]
    end

    plot_nb_eval(rtol_vec, hai, hc) |> display
    plot_error(rtol_vec, hai, hc, Iᵣ) |> display

    return nothing
end
