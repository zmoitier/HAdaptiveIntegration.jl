using CairoMakie
using HAdaptiveIntegration.Rule: GenzMalik, embedded_cubature, orders
using HAdaptiveIntegration: Rectangle, allocate_buffer, integrate
using HCubature
using StaticArrays

include("utils.jl")

const QRULE = GenzMalik{2}()
const DOMAIN = Rectangle((0, 0), (1, 1))
const EC = embedded_cubature(QRULE)
const RTOL_VALUES = [1 / 10^i for i in 1:10]

function add_mesh_inset!(fig_pos, integrand)
    buffer = allocate_buffer(integrand, DOMAIN)
    integrate(integrand, DOMAIN; rtol=1e-6, rule=EC, buffer=buffer)

    ax = Axis(
        fig_pos;
        width=Relative(0.42),
        height=Relative(0.42),
        halign=1.0,
        valign=1.0,
        backgroundcolor=:white,
        aspect=DataAspect(),
    )
    hidedecorations!(ax)
    hidespines!(ax)

    n = 256
    xs = range(0, 1, length=n)
    ys = range(0, 1, length=n)
    z = [integrand(SVector(xi, yi)) for xi in xs, yi in ys]
    heatmap!(ax, xs, ys, z; colormap=:viridis, alpha=0.6, rasterize=16)

    for el in buffer.valtree
        l, h = el[1].corners
        poly!(
            ax,
            [
                Point2f(l[1], l[2]),
                Point2f(h[1], l[2]),
                Point2f(h[1], h[2]),
                Point2f(l[1], h[2])
            ];
            strokecolor=(:white, 0.7),
            strokewidth=0.15,
            color=(:white, 0),
        )
    end
    poly!(
        ax,
        [Point2f(0, 0), Point2f(1, 0), Point2f(1, 1), Point2f(0, 1)];
        strokecolor=:black,
        strokewidth=0.7,
        color=(:black, 0),
    )

    return nothing
end

function reference(integrand)
    I, E = integrate(integrand, DOMAIN; rule=EC, rtol=REFTOL, maxsubdiv=typemax(Int))
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
        I=zeros(length(RTOL_VALUES)),
        E=zeros(length(RTOL_VALUES)),
        N=zeros(length(RTOL_VALUES)),
    )
    hc = (
        I=zeros(length(RTOL_VALUES)),
        E=zeros(length(RTOL_VALUES)),
        N=zeros(length(RTOL_VALUES)),
    )

    for (i, rtol) in enumerate(RTOL_VALUES)
        counter_hai[] = 0
        I, E = integrate(hai_integrand, DOMAIN; rule=EC, rtol=rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter_hai[]

        counter_hc[] = 0
        I, E = hcubature(hc_integrand, zeros(2), ones(2); rtol=rtol)
        hc.I[i], hc.E[i], hc.N[i] = I, E, counter_hc[]
    end

    return hai, hc, Iref, high, low
end

function plot_hai_hc_panel!(ax, hai_data, hc_data, Iref, high, low, order_dim, c_hai, c_hc)
    p1 = scatterlines!(
        ax, hai_data.N, abs.(hai_data.I .- Iref) ./ abs(Iref); color=c_hai, marker=:circle
    )
    p2 = scatterlines!(
        ax, hai_data.N, hai_data.E ./ abs(Iref); color=c_hai, marker=:rect, linestyle=:dash
    )
    p3 = scatterlines!(
        ax, hc_data.N, abs.(hc_data.I .- Iref) ./ abs(Iref); color=c_hc, marker=:circle
    )
    p4 = scatterlines!(
        ax, hc_data.N, hc_data.E ./ abs(Iref); color=c_hc, marker=:rect, linestyle=:dash
    )

    idx_ref = length(hai_data.N)
    N_ref = hai_data.N[idx_ref]
    e_ref = abs(hai_data.I[idx_ref] - Iref) / abs(Iref)
    est_ref = hai_data.E[idx_ref] / abs(Iref)
    p5 = lines!(
        ax,
        hai_data.N,
        e_ref .* (hai_data.N ./ N_ref) .^ (-high / order_dim);
        color=:black,
        linestyle=:dot,
        linewidth=1,
    )
    p6 = lines!(
        ax,
        hai_data.N,
        est_ref .* (hai_data.N ./ N_ref) .^ (-low / order_dim);
        color=:gray,
        linestyle=:dot,
        linewidth=1,
    )
    return p1, p2, p3, p4, p5, p6
end

function main()
    fct_point, fct_sphere, fct_plane = make_features(2)

    println("Running convergence for point feature...")
    hai_point, hc_point, Iref_point, high, low = run_convergence(fct_point)
    println("Running convergence for hypersphere feature...")
    hai_sphere, hc_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
    println("Running convergence for hyperplane feature...")
    hai_plane, hc_plane, Iref_plane, _, _ = run_convergence(fct_plane)

    fig_w, fig_h, font_size = figure_sizes(1, 0.4)
    set_theme!(paper_theme(font_size))
    fig_cvg = Figure(size=(fig_w, fig_h))

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
            axes,
            Axis(
                fig_cvg[1, col];
                xlabel=L"N",
                ylabel=ylab,
                xscale=log10,
                yscale=log10,
                title=title,
            )
        )[end]

        p1, p2, p3, p4, p5, p6 = plot_hai_hc_panel!(
            ax, hai_data, hc_data, Iref, high, low, 2, c_hai, c_hc
        )

        add_mesh_inset!(fig_cvg[1, col], integrand)

        col == 1 && append!(legend_entries, [p5, p6])
    end

    linkaxes!(axes...)
    xlims!(axes[1], 7e1, 8e5)
    ylims!(axes[1], 1e-15, 1e0)

    Legend(
        fig_cvg[2, :],
        [
            Makie.LineElement(color=c_hai),
            Makie.LineElement(color=c_hc),
            Makie.LineElement(color=:black, linestyle=:solid),
            Makie.LineElement(color=:black, linestyle=:dash),
            legend_entries...,
        ],
        [
            "HAI",
            "HC",
            "Actual error",
            "Esti. error",
            L"\mathcal{O}(N^{-%$(high)/2})",
            L"\mathcal{O}(N^{-%$(low)/2})",
        ];
        orientation=:horizontal,
        framevisible=false,
        halign=1,
    )
    rowgap!(fig_cvg.layout, 1, 0)

    save("./paper/image/cvg_rectangle.pdf", fig_cvg; pt_per_unit=1)
    println("Saved cvg_rectangle.pdf")

    return nothing
end

main()
