using GLMakie
using HAdaptiveIntegration.Domain:
    AbstractDomain,
    Cuboid,
    Orthotope,
    Rectangle,
    Segment,
    Simplex,
    Tetrahedron,
    Triangle,
    map_from_reference
using HAdaptiveIntegration.Rule:
    AbstractRule,
    CUBE_BE115,
    CUBE_BE65,
    EmbeddedCubature,
    GenzMalik,
    GrundmannMoeller,
    RadonLaurie,
    SEGMENT_GK15,
    SEGMENT_GK31,
    SEGMENT_GK7,
    SQUARE_CH21,
    SQUARE_CH25,
    TabulatedEmbeddedCubature,
    embedded_cubature,
    orders
using StaticArrays

function vertices(domain::Triangle)
    return Vector(domain.vertices)
end

function vertices(domain::Rectangle)
    a, b = domain.corners
    return [a, SVector(b[1], a[2]), b, SVector(a[1], b[2])]
end

function vertices(domain::Tetrahedron)
    return Vector(domain.vertices)
end

function vertices(domain::Cuboid)
    a, b = domain.corners
    return [
        a,
        SVector(b[1], a[2], a[3]),
        SVector(a[1], b[2], a[3]),
        SVector(b[1], b[2], a[3]),
        SVector(a[1], a[2], b[3]),
        SVector(b[1], a[2], b[3]),
        SVector(a[1], b[2], b[3]),
        b,
    ]
end

edges(::Tetrahedron) = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))

function edges(::Cuboid)
    return (
        (1, 2),
        (1, 3),
        (2, 4),
        (3, 4),
        (5, 6),
        (5, 7),
        (6, 8),
        (7, 8),
        (1, 5),
        (2, 6),
        (3, 7),
        (4, 8),
    )
end

function color_weighs(weights::Vector{T}) where {T<:Real}
    tab10 = Makie.to_colormap(:tab10)
    return [w > 0 ? (tab10[4], 0.75) : (tab10[1], 0.75) for w in weights]
end

function size_weights(weights::Vector{T}, w_max::T) where {T<:Real}
    return 32 * sqrt.(abs.(weights) ./ w_max)
end

function _plot_segment(ec::EmbeddedCubature, dh::Int, dl::Int)
    H, L = length(ec.weights_high), length(ec.weights_low)

    nodes = reinterpret(Float32, ec.nodes)

    fig = Figure()
    ax = Axis(fig[1, 1])
    tab10 = Makie.to_colormap(:tab10)

    scatter!(ax, nodes[1:L], zeros(L); color=tab10[1])
    scatter!(
        ax,
        nodes[1:L],
        ec.weights_low ./ maximum(abs.(ec.weights_low));
        marker=:xcross,
        color=tab10[1],
        label="low order = $dl",
    )

    scatter!(ax, nodes[(L + 1):H], zeros(H - L); color=tab10[2])
    scatter!(
        ax,
        nodes,
        ec.weights_high ./ maximum(abs.(ec.weights_high));
        marker=:cross,
        color=tab10[2],
        label="high order = $dh",
    )

    axislegend(ax)

    return fig
end

function _plot_2d(
    ec::EmbeddedCubature{2}, domain::AbstractDomain{2}, order_high::Int, order_low::Int
)
    L = length(ec.weights_low)
    Φ, _ = map_from_reference(domain)
    nodes = Φ.(ec.nodes)

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
        marker=:xcross,
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

function _plot_3d(
    ec::EmbeddedCubature{3}, domain::AbstractDomain{3}, order_high::Int, order_low::Int
)
    L = length(ec.weights_low)
    Φ, _ = map_from_reference(domain)
    nodes = Φ.(ec.nodes)

    fig = Figure()

    ax_args = Dict(:aspect => :data, :xlabel => "x", :ylabel => "y", :zlabel => "z")
    axs = [
        Axis3(fig[1, 1]; title="high order = $order_high", ax_args...),
        Axis3(fig[1, 2]; title="low order = $order_low", ax_args...),
    ]

    pts = vertices(domain)
    for ax in axs
        for (i, j) in edges(domain)
            p, q = pts[i], pts[j]
            lines!(ax, [p[1], q[1]], [p[2], q[2]], [p[3], q[3]]; color=:black, linewidth=1)
        end
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
        marker=:xcross,
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

function plot_rule(
    rule::AbstractRule{DOM}
) where {DOM<:Union{Segment,Simplex{1},Orthotope{1}}}
    ec = embedded_cubature(rule, Float32)
    oh, ol = orders(rule)
    return _plot_segment(ec, oh, ol)
end

function plot_rule(rule::AbstractRule{DOM}) where {DOM<:Union{Triangle,Simplex{2}}}
    ec = embedded_cubature(rule, Float32)
    oh, ol = orders(rule)
    domain = Triangle{Float32}((1, 0), (-0.5, √3 / 2), (-0.5, -√3 / 2))
    return _plot_2d(ec, domain, oh, ol)
end

function plot_rule(rule::AbstractRule{DOM}) where {DOM<:Union{Rectangle,Orthotope{2}}}
    ec = embedded_cubature(rule, Float32)
    oh, ol = orders(rule)
    domain = Rectangle{Float32}((-1, -1), (1, 1))
    return _plot_2d(ec, domain, oh, ol)
end

function plot_rule(rule::AbstractRule{DOM}) where {DOM<:Union{Tetrahedron,Simplex{3}}}
    ec = embedded_cubature(rule, Float32)
    oh, ol = orders(rule)

    s3, s6 = sqrt(3), sqrt(6)
    domain = Tetrahedron{Float32}(
        (0, 0, 0), (1, 0, 0), (1 / 2, s3 / 2, 0), (1 / 2, s3 / 6, s6 / 3)
    )
    return _plot_3d(ec, domain, oh, ol)
end

function plot_rule(rule::AbstractRule{DOM}) where {DOM<:Union{Cuboid,Orthotope{3}}}
    ec = embedded_cubature(rule, Float32)
    oh, ol = orders(rule)
    domain = Cuboid{Float32}((-1, -1, -1), (1, 1, 1))
    return _plot_3d(ec, domain, oh, ol)
end

function main()
    print("""
    Available rules (# default):
      > Segment:
        - GenzMalik{1}()
        - GrundmannMoeller{1}(deg_high, deg_low)
        - SEGMENT_GK7
        # SEGMENT_GK15
        - SEGMENT_GK31
      > Triangle:
        - GrundmannMoeller{2}(deg_high, deg_low)
        # RadonLaurie()
      > Rectangle:
        - GenzMalik{2}()
        - SQUARE_CH21
        # SQUARE_CH25
      > Tetrahedron:
        # GrundmannMoeller{3}(7, 5)
        - GrundmannMoeller{3}(deg_high, deg_low)
      > Cuboid:
        # CUBE_BE65
        - CUBE_BE115
        - GenzMalik{3}()
    """)

    return nothing
end

main()
