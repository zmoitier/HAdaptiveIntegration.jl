using CairoMakie
using HAdaptiveIntegration.Rule: RadonLaurie, embedded_cubature, orders
using HAdaptiveIntegration: Triangle, allocate_buffer, integrate
using StaticArrays

include("utils.jl")

const QRULE = RadonLaurie()
const DOMAIN = Triangle((0, 0), (1, 0), (0, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:10]

function reference(integrand)
    I, E = integrate(integrand, DOMAIN; rule=EC, rtol=REFTOL, maxsubdiv=typemax(Int))
    return I, E
end

function run_convergence(integrand)
    high, low = orders(QRULE) .+ 1
    counter = Ref(0)
    counted_integrand = (x) -> (counter[] += 1; integrand(x))

    Iref, _ = reference(integrand)

    hai = (
        I=zeros(length(RTOL_VALUES)),
        E=zeros(length(RTOL_VALUES)),
        N=zeros(length(RTOL_VALUES)),
    )

    for (i, rtol) in enumerate(RTOL_VALUES)
        counter[] = 0
        I, E = integrate(counted_integrand, DOMAIN; rule=EC, rtol=rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

function add_mesh_inset!(fig_pos, integrand)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; rtol=1e-6, buffer=buffer)

    ax = Axis(
        fig_pos;
        width=Relative(0.42),
        height=Relative(0.42),
        halign=0,
        valign=0,
        backgroundcolor=:white,
        aspect=DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 256
    xs = range(0, 1, length=n)
    ys = range(0, 1, length=n)
    z = [xi + yi < 1 ? integrand(SVector(xi, yi)) : NaN for xi in xs, yi in ys]
    heatmap!(ax, xs, ys, z; colormap=:viridis, alpha=0.6, rasterize=16)

    for el in buffer.valtree
        v1, v2, v3 = el[1].vertices
        poly!(
            ax,
            [Point2f(v1), Point2f(v2), Point2f(v3)];
            strokecolor=(:white, 0.7),
            strokewidth=0.15,
            color=(:white, 0),
        )
    end
    v1, v2, v3 = DOMAIN.vertices
    poly!(
        ax,
        [Point2f(v1), Point2f(v2), Point2f(v3)];
        strokecolor=:black,
        strokewidth=0.7,
        color=(:black, 0),
    )

    return nothing
end

function plot_hai_panel!(ax, hai_data, Iref, high, low, order_dim, color)
    p1 = scatterlines!(
        ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color=color, marker=:circle
    )
    p2 = scatterlines!(
        ax, hai_data.N, hai_data.E ./ abs(Iref); color=color, marker=:rect, linestyle=:dash
    )
    idx_ref = length(hai_data.N)
    p3 = lines!(
        ax,
        hai_data.N,
        (abs(hai_data.I[idx_ref] - Iref) / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-high / order_dim);
        color=:black,
        linestyle=:dot,
        linewidth=1,
    )
    p4 = lines!(
        ax,
        hai_data.N,
        (hai_data.E[idx_ref] / abs(Iref)) .* (hai_data.N ./ hai_data.N[idx_ref]) .^ (-low / order_dim);
        color=:gray,
        linestyle=:dot,
        linewidth=1,
    )
    return p1, p2, p3, p4
end

function main()
    fct_point, fct_sphere, fct_plane = make_features(2)

    println("Running convergence for point feature...")
    hai_point, Iref_point, high, low = run_convergence(fct_point)
    println("Running convergence for hypersphere feature...")
    hai_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
    println("Running convergence for hyperplane feature...")
    hai_plane, Iref_plane, _, _ = run_convergence(fct_plane)

    fig_w, fig_h, font_size = figure_sizes(1, 0.4)
    set_theme!(paper_theme(font_size))
    fig_cvg = Figure(size=(fig_w, fig_h))

    c_hai = Makie.wong_colors()[1]

    axes_cvg = Axis[]
    legend_plots = Any[]
    for (col, hai_data, Iref, integrand, title) in (
        (1, hai_point, Iref_point, fct_point, "Point Feature"),
        (2, hai_sphere, Iref_sphere, fct_sphere, "Hypersphere Feature"),
        (3, hai_plane, Iref_plane, fct_plane, "Hyperplane Feature"),
    )
        ylab = col == 1 ? "Relative error" : ""
        ax = push!(
            axes_cvg, Axis(
            fig_cvg[1, col];
            xlabel=L"N",
            ylabel=ylab,
            xscale=log10,
            yscale=log10,
            title=title,
        )
        )[end]

        p1, p2, p3, p4 = plot_hai_panel!(ax, hai_data, Iref, high, low, 2, c_hai)

        add_mesh_inset!(fig_cvg[1, col], integrand)

        col == 1 && append!(legend_plots, [p1, p2, p3, p4])
    end

    linkaxes!(axes_cvg...)
    xlims!(axes_cvg[1], 1e2, 1e6)
    ylims!(axes_cvg[1], 1e-15, 1e0)

    Legend(
        fig_cvg[2, :],
        legend_plots,
        [
            "Actual error",
            "Estimated error",
            L"\mathcal{O}(N^{-%$(high)/2})",
            L"\mathcal{O}(N^{-%$(low)/2})",
        ];
        orientation=:horizontal,
        framevisible=false,
        halign=1,
    )
    rowgap!(fig_cvg.layout, 1, 0)

    save("./paper/image/cvg_triangle.pdf", fig_cvg; pt_per_unit=1)
    println("Saved cvg_triangle.pdf")

    return nothing
end

main()
