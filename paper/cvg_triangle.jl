using Quadmath
using HAdaptiveIntegration: integrate, Triangle, allocate_buffer
using HAdaptiveIntegration.Rule: RadonLaurie, embedded_cubature, orders
using CairoMakie
using StaticArrays
using LinearAlgebra

include("util.jl")

function reference(f)
    f_128 = (x) -> f(Float128.(x))
    rule = embedded_cubature(RadonLaurie(), Float128)
    dom = Triangle{Float128}((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    I, E = integrate(f_128, dom; rtol = Float128(1.0e-10), rule = rule)
    return I, E
end

function run_convergence(fct)
    high, low = orders(RadonLaurie()) .+ 1
    counter = Ref(0)
    fct_count = (x) -> (counter[] += 1; fct(x))
    dom = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

    rtol_vec = [1 / 10^i for i in 1:10]
    Iref, _ = reference(fct)

    hai = (I = zeros(length(rtol_vec)), E = zeros(length(rtol_vec)), N = zeros(length(rtol_vec)))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, dom; rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

function add_mesh_inset!(fig_pos, fct)
    dom = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    buffer = allocate_buffer(fct, dom)
    integrate(fct, dom; buffer)

    ax = Axis(
        fig_pos;
        width = Relative(0.42),
        height = Relative(0.42),
        halign = 0.08,
        valign = 0.08,
        backgroundcolor = :white,
        aspect = DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    # heatmap of the integrand
    n = 150
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    z = [xi + yi < 1 ? fct(SVector(xi, yi)) : NaN for xi in xs, yi in ys]
    heatmap!(ax, xs, ys, z; colormap = :viridis, alpha = 0.6)

    # adaptive mesh lines in semi-transparent white
    for el in buffer.valtree
        tri = el[1]
        xt, yt = plot_triangle(tri)
        lines!(ax, xt, yt; color = (:white, 0.25), linewidth = 0.5)
    end
    xt, yt = plot_triangle(dom)
    return lines!(ax, xt, yt; color = :black, linewidth = 2)
end

## definitions
ϵ = 0.025
x₀ = SVector(1 / pi, 1 / pi)
r₀ = 0.5

fct_point = (x) -> scaled_mollifier(norm(x - x₀), ϵ, 2)
fct_curve = (x) -> scaled_mollifier(dot(x, x) - r₀^2, ϵ, 1)
fct_line = (x) -> scaled_mollifier(x[1] - 1 / π, ϵ, 1)

println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for curve feature...")
hai_curve, Iref_curve, _, _ = run_convergence(fct_curve)
println("Running convergence for line feature...")
hai_line, Iref_line, _, _ = run_convergence(fct_line)

fig_cvg = Figure(size = (1200, 500))

axes_cvg = Axis[]
legend_plots = Any[]
for (col, hai, Iref, fct, title) in (
        (1, hai_point, Iref_point, fct_point, "Point Feature"),
        (2, hai_curve, Iref_curve, fct_curve, "Curve Feature"),
        (3, hai_line, Iref_line, fct_line, "Line Feature"),
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
        ax, hai.N, (abs(hai.I[idx_ref] - Iref) / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-high / 2);
        color = :black, linestyle = :dash,
    )
    p4 = lines!(
        ax, hai.N, (hai.E[idx_ref] / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-low / 2);
        color = :gray, linestyle = :dash,
    )

    add_mesh_inset!(fig_cvg[1, col], fct)

    col == 1 && append!(legend_plots, [p1, p2, p3, p4])
end

linkaxes!(axes_cvg...)
Legend(
    fig_cvg[2, :],
    legend_plots,
    ["Actual error", "Estimated error", L"\mathcal{O}(N^{-%$(high)/2})", L"\mathcal{O}(N^{-%$(low)/2})"];
    orientation = :horizontal,
    framevisible = false,
)

save(joinpath(@__DIR__, "cvg_triangle.png"), fig_cvg)
println("Saved cvg_triangle.png")
