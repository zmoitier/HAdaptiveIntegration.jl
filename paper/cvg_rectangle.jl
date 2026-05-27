using CairoMakie
using HAdaptiveIntegration.Rule: GenzMalik, SQUARE_CH25, embedded_cubature, orders
using HAdaptiveIntegration: Rectangle, allocate_buffer, integrate
using HCubature
using LinearAlgebra
using StaticArrays

include("util.jl")

## Quadrature rule used for HAdaptiveIntegration (swap to e.g. SQUARE_CH25 to compare)
const QRULE = GenzMalik{2}()
# const QRULE = SQUARE_CH25
const DOMAIN = Rectangle((0, 0), (1, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:10]

function add_mesh_inset!(fig_pos, integrand)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; buffer, rule = EC)

    ax = Axis(
        fig_pos;
        width = Relative(0.42),
        height = Relative(0.42),
        halign = 0.98,
        valign = 1.0,
        backgroundcolor = :white,
        aspect = DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 150
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    z = [integrand(SVector(xi, yi)) for xi in xs, yi in ys]
    heatmap!(ax, xs, ys, z; colormap = :viridis, alpha = 0.6)

    for el in buffer.valtree
        rect = el[1]
        xr, yr = plot_rectangle(rect)
        lines!(ax, xr, yr; color = (:white, 0.25), linewidth = 0.5)
    end
    return lines!(ax, [0, 1, 1, 0, 0], [0, 0, 1, 1, 0]; color = :black, linewidth = 2)
end

function reference(integrand)
    I, E = integrate(integrand, DOMAIN; rule = EC, rtol = REFTOL, maxsubdiv = typemax(Int))
    return I, E
end

function run_convergence(integrand)
    high, low = orders(QRULE) .+ 1
    counter_hai = Ref(0)
    counter_hc = Ref(0)
    hai_integrand = (x) -> (counter_hai[] += 1; integrand(x))
    hc_integrand = (x) -> (counter_hc[] += 1; integrand(x))

    Iref, _ = reference(integrand)

    hai = (I = zeros(length(RTOL_VALUES)), E = zeros(length(RTOL_VALUES)), N = zeros(length(RTOL_VALUES)))
    hc = (I = zeros(length(RTOL_VALUES)), E = zeros(length(RTOL_VALUES)), N = zeros(length(RTOL_VALUES)))

    for (i, rtol) in enumerate(RTOL_VALUES)
        counter_hai[] = 0
        I, E = integrate(hai_integrand, DOMAIN; rule = EC, rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter_hai[]

        counter_hc[] = 0
        I, E = hcubature(hc_integrand, zeros(2), ones(2); rtol = rtol)
        hc.I[i], hc.E[i], hc.N[i] = I, E, counter_hc[]
    end

    return hai, hc, Iref, high, low
end

function plot_hai_hc_panel!(ax, hai_data, hc_data, Iref, high, low, order_dim, c_hai, c_hc)
    p1 = scatterlines!(ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color = c_hai, marker = :circle)
    p2 = scatterlines!(ax, hai_data.N, hai_data.E ./ abs(Iref); color = c_hai, marker = :rect, linestyle = :dash)
    p3 = scatterlines!(ax, hc_data.N, abs.(hc_data.I .- Iref) ./ abs(Iref); color = c_hc, marker = :circle)
    p4 = scatterlines!(ax, hc_data.N, hc_data.E ./ abs(Iref); color = c_hc, marker = :rect, linestyle = :dash)

    idx_ref = length(hai_data.N)
    N_ref = hai_data.N[idx_ref]
    e_ref = abs(hai_data.I[idx_ref] - Iref) / abs(Iref)
    est_ref = hai_data.E[idx_ref] / abs(Iref)
    p5 = lines!(ax, hai_data.N, e_ref .* (hai_data.N ./ N_ref) .^ (-high / order_dim); color = :black, linestyle = :dot, linewidth = 2)
    p6 = lines!(ax, hai_data.N, est_ref .* (hai_data.N ./ N_ref) .^ (-low / order_dim); color = :gray, linestyle = :dot, linewidth = 2)
    return p1, p2, p3, p4, p5, p6
end

fct_point, fct_sphere, fct_plane = make_features(2)

println("Running convergence for point feature...")
hai_point, hc_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for hypersphere feature...")
hai_sphere, hc_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
println("Running convergence for hyperplane feature...")
hai_plane, hc_plane, Iref_plane, _, _ = run_convergence(fct_plane)

fig = Figure(size = (1300, 520))

c_hai = Makie.wong_colors()[1]
c_hc = Makie.wong_colors()[2]

axes = Axis[]
legend_entries = Any[]
for (col, hai_data, hc_data, Iref, integrand, title) in (
        (1, hai_point, hc_point, Iref_point, fct_point, "Point Feature"),
        (2, hai_sphere, hc_sphere, Iref_sphere, fct_sphere, "Hypersphere Feature"),
        (3, hai_plane, hc_plane, Iref_plane, fct_plane, "Hyperplane Feature"),
    )
    ylab = col == 1 ? "Relative error" : ""
    ax = push!(
        axes, Axis(
            fig[1, col];
            xlabel = "N (number of evaluations)",
            ylabel = ylab,
            xscale = log10,
            yscale = log10,
            title = title,
        )
    )[end]

    p1, p2, p3, p4, p5, p6 = plot_hai_hc_panel!(ax, hai_data, hc_data, Iref, high, low, 2, c_hai, c_hc)

    add_mesh_inset!(fig[1, col], integrand)

    col == 1 && append!(legend_entries, [p1, p2, p3, p4, p5, p6])
end

linkaxes!(axes...)
Legend(
    fig[2, :],
    legend_entries,
    [
        "HAI (actual)", "HAI (estimated)",
        "HCubature (actual)", "HCubature (estimated)",
        L"\mathcal{O}(N^{-%$(high)/2})", L"\mathcal{O}(N^{-%$(low)/2})",
    ];
    orientation = :horizontal,
    framevisible = false,
)
save(joinpath(@__DIR__, "cvg_rectangle.png"), fig)
println("Saved cvg_rectangle.png")
