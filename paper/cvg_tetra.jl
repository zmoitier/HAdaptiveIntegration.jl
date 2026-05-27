using Quadmath
using HAdaptiveIntegration: integrate, Tetrahedron
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature, orders
using CairoMakie
using StaticArrays
using LinearAlgebra

include("util.jl")

function reference(f)
    T = Float64
    f_ext = (x) -> f(T.(x))
    rule = embedded_cubature(GrundmannMoeller{3}(7, 5), T)
    dom = Tetrahedron{T}((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    I, E = integrate(f_ext, dom; rtol = T(1.0e-10), rule = rule, maxsubdiv = 10^6)
    return I, E
end

function run_convergence(fct)
    high, low = orders(GrundmannMoeller{3}(7, 5)) .+ 1
    counter = Ref(0)
    fct_count = (x) -> (counter[] += 1; fct(x))
    dom = Tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    rtol_vec = [1 / 10^i for i in 1:8]
    Iref, _ = reference(fct)

    hai = (I = zeros(length(rtol_vec)), E = zeros(length(rtol_vec)), N = zeros(length(rtol_vec)))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, dom; rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

## definitions
ϵ = 0.05
x₀ = SVector(1 / pi, 1 / pi, 1 / pi)
r₀ = 0.5

fct_point = (x) -> scaled_mollifier(norm(x - x₀), ϵ, 3)
fct_surface = (x) -> scaled_mollifier(dot(x, x) - r₀^2, ϵ, 1)
fct_plane = (x) -> scaled_mollifier(x[1] - 1 / π, ϵ, 1)

println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for surface feature...")
hai_curve, Iref_curve, _, _ = run_convergence(fct_surface)
println("Running convergence for plane feature...")
hai_plane, Iref_plane, _, _ = run_convergence(fct_plane)

fig_cvg = Figure(size = (1200, 500))

axes_cvg = Axis[]
legend_plots = Any[]
for (col, hai, Iref, title) in (
        (1, hai_point, Iref_point, "Point Feature"),
        (2, hai_curve, Iref_curve, "Surface Feature"),
        (3, hai_plane, Iref_plane, "Plane Feature"),
    )
    ylab = col == 1 ? "Relative error" : ""
    ax = push!(
        axes_cvg, Axis(
            fig_cvg[1, col];
            xlabel = "N (number of evaluations)",
            ylabel = ylab,
            xscale = log10,
            yscale = log10,
            title = title,
        )
    )[end]

    p1 = scatterlines!(ax, hai.N, abs.(hai.I .- Iref) ./ abs(Iref); marker = :circle)
    p2 = scatterlines!(ax, hai.N, hai.E ./ abs(Iref); marker = :rect)
    idx_ref = length(hai.N)
    p3 = lines!(
        ax, hai.N, (abs(hai.I[idx_ref] - Iref) / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-high / 3);
        color = :black, linestyle = :dash,
    )
    p4 = lines!(
        ax, hai.N, (hai.E[idx_ref] / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-low / 3);
        color = :gray, linestyle = :dash,
    )

    col == 1 && append!(legend_plots, [p1, p2, p3, p4])
end

linkaxes!(axes_cvg...)
Legend(
    fig_cvg[2, :],
    legend_plots,
    ["Actual error", "Estimated error", L"\mathcal{O}(N^{-%$(high)/3})", L"\mathcal{O}(N^{-%$(low)/3})"];
    orientation = :horizontal,
    framevisible = false,
)

save(joinpath(@__DIR__, "cvg_tetrahedron.png"), fig_cvg)
println("Saved cvg_tetrahedron.png")
