using GLMakie
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature, orders
using HAdaptiveIntegration: Tetrahedron, allocate_buffer, integrate
using LinearAlgebra
using StaticArrays

include("util.jl")
# Force a 3D-capable backend to avoid Cairo volume limitations.
GLMakie.activate!()

const QRULE = GrundmannMoeller{3}(7, 5)
const DOMAIN = Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:8]

function plot_tetrahedron_edges(tet)
    v = tet.vertices
    edges = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    for (i, j) in edges
        push!(xs, v[i][1], v[j][1], NaN)
        push!(ys, v[i][2], v[j][2], NaN)
        push!(zs, v[i][3], v[j][3], NaN)
    end
    return xs, ys, zs
end

function add_mesh_inset!(fig_pos, integrand, clip::Bool)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; buffer, rtol = 1.0e-2)

    ax = Axis3(
        fig_pos;
        width = Relative(0.42),
        height = Relative(0.42),
        halign = 0.08,
        valign = 0.08,
        aspect = :data,
        backgroundcolor = :white,
        azimuth = pi / 3,
        elevation = pi / 6,
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 45
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    zs = range(0, 1, length = n)
    vals = if clip
        [xi + yi + zi < 1 ? integrand(SVector(xi, yi, zi)) : NaN for xi in xs, yi in ys, zi in zs]
    else
        [integrand(SVector(xi, yi, zi)) for xi in xs, yi in ys, zi in zs]
    end
    finite_vals = filter(isfinite, vec(vals))
    vmin, vmax = extrema(finite_vals)
    levels = collect(range(vmin + 0.2 * (vmax - vmin), vmin + 0.8 * (vmax - vmin), length = 4))

    for el in buffer.valtree
        tet = el[1]
        xe, ye, ze = plot_tetrahedron_edges(tet)
        # Keep adaptive cell edges visible above translucent contour surfaces.
        lines!(ax, xe, ye, ze; color = (:gray, 0.75), linewidth = 0.75, overdraw = true)
    end

    xe, ye, ze = plot_tetrahedron_edges(DOMAIN)
    lines!(ax, xe, ye, ze; color = :black, linewidth = 2, overdraw = true)

    contour!(
        ax,
        first(xs) .. last(xs),
        first(ys) .. last(ys),
        first(zs) .. last(zs),
        vals;
        levels = levels,
        colormap = :viridis,
        alpha = 0.25,
        overdraw = true,
    )

    return ax
end

function reference(integrand)
    I, E = integrate(integrand, DOMAIN; rule = EC, rtol = REFTOL, maxsubdiv = typemax(Int))
    return I, E
end

function run_convergence(integrand)
    high, low = orders(QRULE) .+ 1
    counter = Ref(0)
    counted_integrand = (x) -> (counter[] += 1; integrand(x))

    Iref, _ = reference(integrand)

    hai = (I = zeros(length(RTOL_VALUES)), E = zeros(length(RTOL_VALUES)), N = zeros(length(RTOL_VALUES)))

    for (i, rtol) in enumerate(RTOL_VALUES)
        counter[] = 0
        I, E = integrate(counted_integrand, DOMAIN; rule = EC, rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

function plot_hai_panel!(ax, hai_data, Iref, high, low, order_dim, color)
    p1 = scatterlines!(ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color = color, marker = :circle)
    p2 = scatterlines!(ax, hai_data.N, hai_data.E ./ abs(Iref); color = color, marker = :rect, linestyle = :dash)
    idx_ref = length(hai_data.N)
    p3 = lines!(
        ax,
        hai_data.N,
        (abs(hai_data.I[idx_ref] - Iref) / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-high / order_dim);
        color = :black,
        linestyle = :dot,
        linewidth = 2,
    )
    p4 = lines!(
        ax,
        hai_data.N,
        (hai_data.E[idx_ref] / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-low / order_dim);
        color = :gray,
        linestyle = :dot,
        linewidth = 2,
    )
    return p1, p2, p3, p4
end

fct_point, fct_sphere, fct_plane = make_features(3)

println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for hypersphere feature...")
hai_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
println("Running convergence for hyperplane feature...")
hai_plane, Iref_plane, _, _ = run_convergence(fct_plane)

fig_cvg = Figure(size = (1200, 500))

c_hai = Makie.wong_colors()[1]

axes_cvg = Axis[]
legend_plots = Any[]
for (col, hai_data, Iref, integrand, title, clip) in (
        (1, hai_point, Iref_point, fct_point, "Point Feature", false),
        (2, hai_sphere, Iref_sphere, fct_sphere, "Hypersphere Feature", false),
        (3, hai_plane, Iref_plane, fct_plane, "Hyperplane Feature", true),
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

    p1, p2, p3, p4 = plot_hai_panel!(ax, hai_data, Iref, high, low, 3, c_hai)

    add_mesh_inset!(fig_cvg[1, col], integrand, clip)

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
