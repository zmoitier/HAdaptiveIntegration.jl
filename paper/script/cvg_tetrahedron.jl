using CairoMakie
using CairoMakie.Makie.GeometryBasics: Mesh, TriangleFace, Point3f
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature, orders
using HAdaptiveIntegration: Tetrahedron, allocate_buffer, integrate
using StaticArrays

include("utils.jl")

const QRULE = GrundmannMoeller{3}(7, 5)
const DOMAIN = Tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:8]

function tet_mesh(tet)
    v = tet.vertices
    pts = [Point3f(v[1]), Point3f(v[2]), Point3f(v[3]), Point3f(v[4])]
    faces = [TriangleFace(1, 2, 3), TriangleFace(1, 2, 4), TriangleFace(1, 3, 4), TriangleFace(2, 3, 4)]
    return Mesh(pts, faces)
end

function add_mesh_inset!(fig_pos, integrand, rtol)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; rtol = rtol, buffer = buffer)

    δ = -0.5
    ax = Axis3(
        fig_pos;
        width = Relative(0.5),
        height = Relative(0.5),
        halign = 0,
        valign = 0,
        aspect = :data,
        backgroundcolor = :white,
        azimuth = pi / 3,
        elevation = pi / 6,
        protrusions = 0,
        limits = (δ, 1, δ, 1, δ, 1),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 64
    xs = range(0, 1, length = n)
    ys = range(0, 1, length = n)
    zs = range(0, 1, length = n)
    vals = [
        xi + yi + zi < 1 ? integrand(SVector(xi, yi, zi)) : NaN
            for xi in xs, yi in ys, zi in zs
    ]

    # Max-intensity projections onto the three coordinate planes (NaN outside simplex).
    vals_finite = replace(vals, NaN => -Inf)
    proj_xy = dropdims(maximum(vals_finite; dims = 3); dims = 3)  # z = 0 plane
    proj_yz = dropdims(maximum(vals_finite; dims = 1); dims = 1)  # x = 0 plane
    proj_xz = dropdims(maximum(vals_finite; dims = 2); dims = 2)  # y = 0 plane
    replace!(proj_xy, -Inf => NaN)
    replace!(proj_yz, -Inf => NaN)
    replace!(proj_xz, -Inf => NaN)

    for el in buffer.valtree
        wireframe!(ax, tet_mesh(el[1]); color = :gray, linewidth = 0.5, overdraw = true)
    end
    wireframe!(ax, tet_mesh(DOMAIN); color = :black, linewidth = 1, overdraw = true)

    # CairoMakie cannot render 3D Volume plots; project the integrand onto the
    # x = 0, y = 0, and z = 0 planes as heatmaps instead.  Offset each plane
    # slightly negative so the three projections sit apart and don't overlap.
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
    counter = Ref(0)
    counted_integrand = (x) -> (counter[] += 1; integrand(x))

    Iref, _ = reference(integrand)

    hai = (
        I = zeros(length(RTOL_VALUES)),
        E = zeros(length(RTOL_VALUES)),
        N = zeros(length(RTOL_VALUES)),
    )

    for (i, rtol) in enumerate(RTOL_VALUES)
        counter[] = 0
        I, E = integrate(counted_integrand, DOMAIN; rule = EC, rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

function plot_hai_panel!(ax, hai_data, Iref, high, low, order_dim, color)
    p1 = scatterlines!(
        ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color = color, marker = :circle
    )
    p2 = scatterlines!(
        ax, hai_data.N, hai_data.E ./ abs(Iref); color = color, marker = :rect, linestyle = :dash
    )
    idx_ref = length(hai_data.N)
    p3 = lines!(
        ax,
        hai_data.N,
        (abs(hai_data.I[idx_ref] - Iref) / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-high / order_dim);
        color = :black,
        linestyle = :dot,
        linewidth = 0.7,
    )
    p4 = lines!(
        ax,
        hai_data.N,
        (hai_data.E[idx_ref] / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-low / order_dim);
        color = :gray,
        linestyle = :dot,
        linewidth = 0.7,
    )
    return p1, p2, p3, p4
end

function main()
    fct_point, fct_sphere, fct_plane = make_features(3)

    println("Running convergence for point feature...")
    hai_point, Iref_point, high, low = run_convergence(fct_point)
    println("Running convergence for hypersphere feature...")
    hai_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
    println("Running convergence for hyperplane feature...")
    hai_plane, Iref_plane, _, _ = run_convergence(fct_plane)

    fig_w, fig_h, font_size = figure_sizes(1, 0.4)
    set_theme!(paper_theme(font_size))
    fig_cvg = Figure(size = (fig_w, fig_h))

    c_hai = Makie.wong_colors()[1]

    axes_cvg = Axis[]
    legend_plots = Any[]
    for (col, hai_data, Iref, integrand, title, rtol) in (
            (1, hai_point, Iref_point, fct_point, "Point Feature", 1.0e-1),
            (2, hai_sphere, Iref_sphere, fct_sphere, "Hypersphere Feature", 4.0e-2),
            (3, hai_plane, Iref_plane, fct_plane, "Hyperplane Feature", 2.0e-2),
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

        p1, p2, p3, p4 = plot_hai_panel!(ax, hai_data, Iref, high, low, 3, c_hai)

        add_mesh_inset!(fig_cvg[1, col], integrand, rtol)

        col == 1 && append!(legend_plots, [p1, p2, p3, p4])
    end

    linkaxes!(axes_cvg...)
    xlims!(axes_cvg[1], 2.0e1, 4.0e7)
    ylims!(axes_cvg[1], 3.0e-12, 4.0e0)

    Legend(
        fig_cvg[2, :],
        legend_plots,
        [
            "Actual error",
            "Estimated error",
            L"\mathcal{O}(N^{-%$(high)/3})",
            L"\mathcal{O}(N^{-%$(low)/3})",
        ];
        orientation = :horizontal,
        framevisible = false,
        halign = 1,
    )
    rowgap!(fig_cvg.layout, 1, 0)

    save("./paper/image/cvg_tetrahedron.pdf", fig_cvg; pt_per_unit = 1)
    println("Saved cvg_tetrahedron.pdf")

    return nothing
end

main()
