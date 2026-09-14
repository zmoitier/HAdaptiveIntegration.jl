using GLMakie
using HAdaptiveIntegration.Rule: CUBE_BE65, GenzMalik, embedded_cubature, orders
using HAdaptiveIntegration: Cuboid, allocate_buffer, integrate
using HCubature
using LinearAlgebra
using StaticArrays

include("util.jl")
# Force a 3D-capable backend to avoid Cairo volume limitations.
GLMakie.activate!()

## Quadrature rule used for HAdaptiveIntegration (swap to GenzMalik to compare)
const QRULE = CUBE_BE65
# const QRULE = GenzMalik{3}()
const DOMAIN = Cuboid((0, 0, 0), (1, 1, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:8]

function plot_cuboid_edges(cub)
    l, h = cub.corners[1], cub.corners[2]
    verts = (
        (l[1], l[2], l[3]),
        (h[1], l[2], l[3]),
        (h[1], h[2], l[3]),
        (l[1], h[2], l[3]),
        (l[1], l[2], h[3]),
        (h[1], l[2], h[3]),
        (h[1], h[2], h[3]),
        (l[1], h[2], h[3]),
    )
    edges = (
        (1, 2), (2, 3), (3, 4), (4, 1),
        (5, 6), (6, 7), (7, 8), (8, 5),
        (1, 5), (2, 6), (3, 7), (4, 8),
    )

    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    for (i, j) in edges
        push!(xs, verts[i][1], verts[j][1], NaN)
        push!(ys, verts[i][2], verts[j][2], NaN)
        push!(zs, verts[i][3], verts[j][3], NaN)
    end
    return xs, ys, zs
end

function add_mesh_inset!(fig_pos, integrand)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; buffer, rtol = 1.0e-2)

    ax = Axis3(
        fig_pos;
        width = Relative(0.42),
        height = Relative(0.42),
        halign = 0.98,
        valign = 1.0,
        aspect = :data,
        backgroundcolor = :white,
        azimuth = -pi * 5 / 12,
        elevation = pi / 6,
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 45
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    zs = range(0, 1, length = n)
    vals = [integrand(SVector(x, y, z)) for x in xs, y in ys, z in zs]
    vmin, vmax = extrema(vals)
    levels = collect(range(vmin + 0.2 * (vmax - vmin), vmin + 0.8 * (vmax - vmin), length = 4))

    for el in buffer.valtree
        cub = el[1]
        xe, ye, ze = plot_cuboid_edges(cub)
        # Keep adaptive cell edges visible above translucent contour surfaces.
        lines!(ax, xe, ye, ze; color = (:gray, 0.5), linewidth = 0.75, overdraw = true)
    end

    xe, ye, ze = plot_cuboid_edges(DOMAIN)
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
        overdraw = true
    )

    return ax
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
        I, E = hcubature(hc_integrand, zeros(3), ones(3); rtol = rtol)
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

fct_point, fct_sphere, fct_plane = make_features(3)

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

    p1, p2, p3, p4, p5, p6 = plot_hai_hc_panel!(ax, hai_data, hc_data, Iref, high, low, 3, c_hai, c_hc)

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
        L"\mathcal{O}(N^{-%$(high)/3})", L"\mathcal{O}(N^{-%$(low)/3})",
    ];
    orientation = :horizontal,
    framevisible = false,
)

outpath = joinpath(@__DIR__, "cvg_cube.png")
save(outpath, fig)
println("Saved $(outpath)")
