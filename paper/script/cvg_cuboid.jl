using CairoMakie
using HAdaptiveIntegration.Rule: CUBE_BE65, embedded_cubature, orders
using HAdaptiveIntegration: Cuboid, allocate_buffer, integrate
using HCubature
using StaticArrays

include("utils.jl")

const QRULE = CUBE_BE65
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

function add_mesh_inset!(fig_pos, integrand, rtol)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; rtol = rtol, buffer = buffer)

    δ = -0.5
    ax = Axis3(
        fig_pos;
        width = Relative(0.5),
        height = Relative(0.5),
        halign = 1,
        valign = 1,
        aspect = :data,
        backgroundcolor = :white,
        azimuth = pi / 3,
        elevation = pi / 6,
        protrusions = 0,
        limits = (δ, 1, δ, 1, δ, 1),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 45
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    zs = range(0, 1, length = n)
    vals = [integrand(SVector(x, y, z)) for x in xs, y in ys, z in zs]

    # Max-intensity projections onto the three coordinate planes.
    proj_xy = dropdims(maximum(vals; dims = 3); dims = 3)  # z = 0 plane
    proj_yz = dropdims(maximum(vals; dims = 1); dims = 1)  # x = 0 plane
    proj_xz = dropdims(maximum(vals; dims = 2); dims = 2)  # y = 0 plane

    for el in buffer.valtree
        cub = el[1]
        xe, ye, ze = plot_cuboid_edges(cub)
        lines!(ax, xe, ye, ze; color = :gray, linewidth = 0.5, overdraw = true)
    end
    xe, ye, ze = plot_cuboid_edges(DOMAIN)
    lines!(ax, xe, ye, ze; color = :black, linewidth = 1, overdraw = true)

    # CairoMakie cannot render 3D Volume plots; project the integrand onto the
    # x = 0, y = 0, and z = 0 planes as heatmaps instead.
    heatmap!(
        ax,
        xs,
        ys,
        proj_xy;
        colormap = :viridis,
        alpha = 0.5,
        overdraw = true,
        transformation = (:xy, δ)
    )
    heatmap!(
        ax,
        ys,
        zs,
        proj_yz;
        colormap = :viridis,
        alpha = 0.5,
        overdraw = true,
        transformation = (:yz, δ)
    )
    heatmap!(
        ax,
        xs,
        zs,
        proj_xz;
        colormap = :viridis,
        alpha = 0.5,
        overdraw = true,
        transformation = (:xz, δ)
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

    hai = (
        I = zeros(length(RTOL_VALUES)),
        E = zeros(length(RTOL_VALUES)),
        N = zeros(length(RTOL_VALUES)),
    )
    hc = (
        I = zeros(length(RTOL_VALUES)),
        E = zeros(length(RTOL_VALUES)),
        N = zeros(length(RTOL_VALUES)),
    )

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
    p1 = scatterlines!(
        ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color = c_hai, marker = :circle
    )
    p2 = scatterlines!(
        ax, hai_data.N, hai_data.E ./ abs(Iref); color = c_hai, marker = :rect, linestyle = :dash
    )
    p3 = scatterlines!(
        ax, hc_data.N, abs.(hc_data.I .- Iref) ./ abs(Iref); color = c_hc, marker = :circle
    )
    p4 = scatterlines!(
        ax, hc_data.N, hc_data.E ./ abs(Iref); color = c_hc, marker = :rect, linestyle = :dash
    )

    idx_ref = length(hai_data.N)
    N_ref = hai_data.N[idx_ref]
    e_ref = abs(hai_data.I[idx_ref] - Iref) / abs(Iref)
    est_ref = hai_data.E[idx_ref] / abs(Iref)
    p5 = lines!(
        ax,
        hai_data.N,
        e_ref .* (hai_data.N ./ N_ref) .^ (-high / order_dim);
        color = :black,
        linestyle = :dot,
        linewidth = 0.7,
    )
    p6 = lines!(
        ax,
        hai_data.N,
        est_ref .* (hai_data.N ./ N_ref) .^ (-low / order_dim);
        color = :gray,
        linestyle = :dot,
        linewidth = 0.7,
    )

    return p1, p2, p3, p4, p5, p6
end

function main()
    fct_point, fct_sphere, fct_plane = make_features(3)

    println("Running convergence for point feature...")
    hai_point, hc_point, Iref_point, high, low = run_convergence(fct_point)
    println("Running convergence for hypersphere feature...")
    hai_sphere, hc_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
    println("Running convergence for hyperplane feature...")
    hai_plane, hc_plane, Iref_plane, _, _ = run_convergence(fct_plane)

    fig_w, fig_h, font_size = figure_sizes(1, 0.4)
    set_theme!(paper_theme(font_size))
    fig_cvg = Figure(size = (fig_w, fig_h))

    c_hai = Makie.wong_colors()[1]
    c_hc = Makie.wong_colors()[2]

    axes_cvg = Axis[]
    legend_plots = Any[]
    for (col, hai_data, hc_data, Iref, integrand, title, rtol) in (
            (1, hai_point, hc_point, Iref_point, fct_point, "Point Feature", 1.0e-1),
            (2, hai_sphere, hc_sphere, Iref_sphere, fct_sphere, "Hypersphere Feature", 3.2e-2),
            (3, hai_plane, hc_plane, Iref_plane, fct_plane, "Hyperplane Feature", 1.0e-2),
        )
        ylab = col == 1 ? "Relative error" : ""
        ax = push!(
            axes_cvg, Axis(
                fig_cvg[1, col];
                xlabel = L"N",
                ylabel = ylab,
                xscale = log10,
                yscale = log10,
                title = title,
            )
        )[end]

        p1, p2, p3, p4, p5, p6 = plot_hai_hc_panel!(
            ax, hai_data, hc_data, Iref, high, low, 3, c_hai, c_hc
        )

        add_mesh_inset!(fig_cvg[1, col], integrand, rtol)

        col == 1 && append!(legend_plots, [p5, p6])
    end

    linkaxes!(axes_cvg...)
    xlims!(axes_cvg[1], 1.0e2, 2.0e7)
    ylims!(axes_cvg[1], 3.0e-12, 1.0e0)

    Legend(
        fig_cvg[2, :],
        [
            Makie.LineElement(color = c_hai),
            Makie.LineElement(color = c_hc),
            Makie.LineElement(color = :black, linestyle = :solid),
            Makie.LineElement(color = :black, linestyle = :dash),
            legend_plots...,
        ],
        [
            "HAI",
            "HC",
            "Actual error",
            "Esti. error",
            L"\mathcal{O}(N^{-%$(high)/3})",
            L"\mathcal{O}(N^{-%$(low)/3})",
        ];
        orientation = :horizontal,
        framevisible = false,
        halign = 1,
    )
    rowgap!(fig_cvg.layout, 1, 0)

    save("./paper/image/cvg_cuboid.pdf", fig_cvg; pt_per_unit = 1)
    println("Saved cvg_cuboid.pdf")

    return nothing
end

main()
