using GLMakie
using HAdaptiveIntegration:
    AbstractDomain,
    AbstractRule,
    EmbeddedCubature,
    Rectangle,
    SEGMENT_GK15,
    SEGMENT_GK31,
    SQUARE_CHG21,
    SQUARE_CHG25,
    Segment,
    TRIANGLE_GM20,
    TRIANGLE_RL19,
    TabulatedEmbeddedCubature,
    Triangle,
    embedded_cubature,
    map_from_reference,
    rectangle,
    triangle
using StaticArrays

function degree(tec::TabulatedEmbeddedCubature)
    return (tec.order_high, tec.order_low)
end

function color_weighs(weights::Vector{Float64})
    tab10 = Makie.to_colormap(:tab10)
    return [w > 0 ? (tab10[4], 0.75) : (tab10[1], 0.75) for w in weights]
end

function size_weights(weights::Vector{Float64}, w_max::Float64)
    return 32 * sqrt.(abs.(weights) ./ w_max)
end

function plot_rule(rule::AbstractRule{Segment})
    ec = embedded_cubature(rule)
    H, L = length(ec.weights_high), length(ec.weights_low)
    dh, dl = degree(rule)

    nodes = reinterpret(Float64, ec.nodes)

    fig = Figure()
    ax = Axis(fig[1, 1]; title=rule.description)
    tab10 = Makie.to_colormap(:tab10)

    scatter!(ax, nodes[1:L], zeros(L); color=tab10[1])
    scatter!(
        ax,
        nodes[1:L],
        ec.weights_low ./ maximum(ec.weights_low);
        marker=:xcross,
        color=tab10[1],
        label="low order = $dl",
    )

    scatter!(ax, nodes[(L + 1):H], zeros(H - L); color=tab10[2])
    scatter!(
        ax,
        nodes,
        ec.weights_high ./ maximum(ec.weights_high);
        marker=:cross,
        color=tab10[2],
        label="high order = $dh",
    )

    axislegend(ax)

    return fig
end

function vertices(domain::Triangle)
    return Vector(domain.vertices)
end

function vertices(domain::Rectangle)
    return [
        domain.low_corner,
        SVector(domain.high_corner[1], domain.low_corner[2]),
        domain.high_corner,
        SVector(domain.low_corner[1], domain.high_corner[2]),
    ]
end

function plot(
    ec::EmbeddedCubature{2}, domain::AbstractDomain{2}, order_high::Int, order_low::Int
)
    L = length(ec.weights_low)
    Φ = map_from_reference(domain)
    nodes = Φ.(ec.nodes)

    display(nodes[1:L])
    display(nodes[(L + 1):end])

    fig = Figure()

    ax_args = Dict(:aspect => 1, :xlabel => "x", :ylabel => "y")
    axs = [
        Axis(fig[1, 1]; title="high order = $order_high", ax_args...),
        Axis(fig[1, 2]; title="low order = $order_low", ax_args...),
    ]

    for ax in axs
        poly!(ax, vertices(domain); :color => (:black, 0.05))
    end

    W = max(maximum(abs.(ec.weights_high)), maximum(abs.(ec.weights_low)))
    scatter!(
        axs[1],
        nodes[1:L];
        color=color_weighs(ec.weights_high[1:L]),
        markersize=size_weights(ec.weights_high[1:L], W),
    )
    scatter!(
        axs[1],
        nodes[(L + 1):end];
        color=color_weighs(ec.weights_high[(L + 1):end]),
        marker=:diamond,
        markersize=size_weights(ec.weights_high[(L + 1):end], W),
    )
    scatter!(
        axs[2],
        nodes[1:L];
        color=color_weighs(ec.weights_low),
        markersize=size_weights(ec.weights_low, W),
    )

    return fig
end

function plot_rule(rule::AbstractRule{Triangle})
    ec = embedded_cubature(rule)
    dh, dl = degree(rule)
    domain = triangle((1, 0), (-0.5, √3 / 2), (-0.5, -√3 / 2))

    return plot(ec, domain, dh, dl)
end

function plot_rule(rule::AbstractRule{Rectangle})
    ec = embedded_cubature(rule)
    dh, dl = degree(rule)
    domain = rectangle((-1, -1), (1, 1))

    return plot(ec, domain, dh, dl)
end

function main()
    # fig = plot_rule(SEGMENT_GK15)
    fig = plot_rule(SQUARE_CHG21)
    # fig = plot_rule(TRIANGLE_GM20)

    display(fig)

    return nothing
end

main()
