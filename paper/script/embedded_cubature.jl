using CairoMakie
using HAdaptiveIntegration: Rectangle, Triangle, default_rule
using HAdaptiveIntegration.Domain: reference_domain

const CMAP = reverse(cgrad(:managua))

function weight_colors(weights, wmin, wmax)
    return [CMAP[(w - wmin) / (wmax - wmin)] for w in weights]
end

function domain_outline(tri::Triangle)
    a, b, c = tri.vertices
    return [Point2f(a[1], a[2]), Point2f(b[1], b[2]), Point2f(c[1], c[2]), Point2f(a[1], a[2])]
end

function domain_outline(rect::Rectangle)
    a, b = rect.corners
    return [Point2f(a[1], a[2]), Point2f(b[1], a[2]), Point2f(b[1], b[2]), Point2f(a[1], b[2]), Point2f(a[1], a[2])]
end

function node_points(nodes)
    return [Point2f(n[1], n[2]) for n in nodes]
end

function plot_cubature!(fig, col, DOM)
    dom = reference_domain(DOM)
    ec = default_rule(dom)
    L = length(ec.weights_low)

    wmax = max(maximum(abs, ec.weights_high), maximum(abs, ec.weights_low))
    wmin = -wmax

    ax = Axis(fig[1, col]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)

    # Domain outline.
    lines!(ax, domain_outline(dom); color = :black, linewidth = 1)

    # High order cubature: x-cross markers.
    scatter!(
        ax, node_points(ec.nodes);
        marker = :xcross,
        color = weight_colors(ec.weights_high, wmin, wmax),
        markersize = 12,
    )

    # Low order cubature: circle markers (reuses the first L nodes).
    scatter!(
        ax, node_points(ec.nodes[1:L]);
        marker = :circle,
        color = :transparent,
        strokecolor = weight_colors(ec.weights_low, wmin, wmax),
        strokewidth = 2,
        markersize = 20,
    )

    Colorbar(
        fig[1, col + 1];
        colormap = CMAP,
        colorrange = (wmin, wmax),
    )

    return wmin, wmax
end

function main()
    fig = Figure(size = (800, 400))

    plot_cubature!(fig, 1, Triangle)
    plot_cubature!(fig, 3, Rectangle)

    Legend(
        fig[2, :],
        [
            [Makie.MarkerElement(marker = :xcross, markersize = 12)],
            [Makie.MarkerElement(marker = :circle, color = :transparent, strokecolor = :black, strokewidth = 2, markersize = 20)],
        ],
        ["High order", "Low order"];
        labelsize = 16,
        orientation = :horizontal,
    )

    save("./paper/image/embedded_cubature.pdf", fig)

    return nothing
end

main()
