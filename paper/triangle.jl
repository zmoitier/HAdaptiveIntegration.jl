using Quadmath
using HAdaptiveIntegration: integrate, Triangle, allocate_buffer
using HAdaptiveIntegration.Rule: RadonLaurie, embedded_cubature, orders
using CairoMakie
using StaticArrays
using LinearAlgebra

include("util.jl")

function plot_mesh_to_axis!(ax, fct, title_str)
    dom = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    buffer = allocate_buffer(fct, dom)
    I, E = integrate(fct, dom; buffer)

    ## plot heatmap
    n = 200
    x = range(0, 1, length = n)
    y = range(0, 1, length = n)
    z = [xi + yi < 1 ? fct(SVector(xi, yi)) : NaN for xi in x, yi in y]
    color = Makie.wong_colors()[2]
    hm = heatmap!(ax, x, y, z; colormap = :viridis)

    ## plot the adaptive mesh
    for el in buffer.valtree
        tri = el[1]
        xt, yt = plot_triangle(tri)
        lines!(ax, xt, yt; color = color, linewidth = 0.5, alpha = 0.3)
    end

    ## plot the original domain
    xt, yt = plot_triangle(dom)
    lines!(ax, xt, yt; color = :black, linewidth = 3)
    ax.title = title_str
    return hm
end

function reference(f)
    # Use Float128 for speed and precision
    f_128 = (x) -> f(Float128.(x))
    rule = embedded_cubature(RadonLaurie(), Float128)
    dom = Triangle{Float128}((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))
    # Tighter tolerance and high maxsubdiv to get a solid ground truth
    I, E = integrate(f_128, dom; rtol = Float128(1.0e-10), rule = rule)
    return I, E
end

function run_convergence(fct)
    high, low = orders(RadonLaurie()) .+ 1
    counter = Ref(0)
    fct_count = (x) -> (counter[] += 1; fct(x))
    dom = Triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

    rtol_vec = [1 / 10^i for i in 1:10]
    Iref, Eref = reference(fct)

    hai = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, dom; rtol = rtol)
        hai.I[i] = I
        hai.E[i] = E
        hai.N[i] = counter[]
    end

    return hai, Iref, high, low
end

## definitions
ϵ = 0.025
x₀ = SVector(1 / 3, 1 / 5)
r₀ = 0.5

# Standardized integrands for mesh plot and convergence
fct_point = (x) -> scaled_mollifier(norm(x - x₀), ϵ, 2)
fct_curve = (x) -> scaled_mollifier(abs(norm(x) - r₀), ϵ, 2)

# 1. Generate Mesh Plots
println("Generating triangle meshes...")
fig_mesh = Figure(size = (800, 350))

ax_mesh1 = Axis(fig_mesh[1, 1], xlabel = "x", ylabel = "y", aspect = DataAspect())
hm1 = plot_mesh_to_axis!(ax_mesh1, fct_point, "Point Feature")

ax_mesh2 = Axis(fig_mesh[1, 2], xlabel = "x", ylabel = "y", aspect = DataAspect())
hm2 = plot_mesh_to_axis!(ax_mesh2, fct_curve, "Curve Feature")

Colorbar(fig_mesh[1, 3], hm1; label = "f(x)", width = 15)

save("triangle_meshes.png", fig_mesh)
println("Saved triangle_meshes.png")

# 2. Generate Convergence Plots
println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for curve feature...")
hai_curve, Iref_curve, _, _ = run_convergence(fct_curve)

fig_cvg = Figure(size = (800, 400))

# Point Subplot
ax_cvg1 = Axis(
    fig_cvg[1, 1], xlabel = "N (number of evaluations)", ylabel = "Relative error",
    xscale = log10, yscale = log10, title = "Point Feature"
)
scatterlines!(ax_cvg1, hai_point.N, abs.(hai_point.I .- Iref_point) ./ abs(Iref_point); marker = :circle, label = "Actual error")
scatterlines!(ax_cvg1, hai_point.N, hai_point.E ./ abs(Iref_point); marker = :rect, label = "Estimated error")

# dashed lines for slopes
lines!(
    ax_cvg1, hai_point.N, (abs(hai_point.I[end] - Iref_point) / abs(Iref_point)) .* (hai_point.N ./ hai_point.N[end]) .^ (-high / 2);
    color = :black, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(high)/2})"
)
lines!(
    ax_cvg1, hai_point.N, (hai_point.E[end] / abs(Iref_point)) .* (hai_point.N ./ hai_point.N[end]) .^ (-low / 2);
    color = :gray, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(low)/2})"
)

axislegend(ax_cvg1; position = :lb)

# Curve Subplot
ax_cvg2 = Axis(
    fig_cvg[1, 2], xlabel = "N (number of evaluations)",
    xscale = log10, yscale = log10, title = "Curve Feature"
)
scatterlines!(ax_cvg2, hai_curve.N, abs.(hai_curve.I .- Iref_curve) ./ abs(Iref_curve); marker = :circle, label = "Actual error")
scatterlines!(ax_cvg2, hai_curve.N, hai_curve.E ./ abs(Iref_curve); marker = :rect, label = "Estimated error")

lines!(
    ax_cvg2, hai_curve.N, (abs(hai_curve.I[end] - Iref_curve) / abs(Iref_curve)) .* (hai_curve.N ./ hai_curve.N[end]) .^ (-high / 2);
    color = :black, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(high)/2})"
)
lines!(
    ax_cvg2, hai_curve.N, (hai_curve.E[end] / abs(Iref_curve)) .* (hai_curve.N ./ hai_curve.N[end]) .^ (-low / 2);
    color = :gray, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(low)/2})"
)

axislegend(ax_cvg2; position = :lb)

save("cvg_triangle.png", fig_cvg)
println("Saved cvg_triangle.png")
