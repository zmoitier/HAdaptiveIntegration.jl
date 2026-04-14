using BenchmarkTools
using GLMakie
using HAdaptiveIntegration: HAdaptiveIntegration as HAI
using HCubature: HCubature as HC
using Printf: @sprintf

function _compute_buffers_and_stats(f, a, b, rtol::Real)
    hai_domain = HAI.Orthotope(a, b)
    hai_buffer = HAI.allocate_buffer(f, hai_domain)
    hai_I, hai_E = HAI.integrate(
        f, hai_domain; buffer=hai_buffer, rtol=rtol, maxsubdiv=typemax(Int)
    )
    hai_E_range = extrema(E for (_, _, E) in hai_buffer.valtree)

    hc_buffer = HC.hcubature_buffer(f, a, b)
    hc_I, hc_E = HC.hcubature(f, a, b; buffer=hc_buffer, rtol=rtol)
    hc_E_range = extrema(box.E for box in hc_buffer.valtree)

    diff_abs = abs(hai_I - hc_I)
    diff_rel = diff_abs / abs((hc_I + hai_I) / 2)
    lims_err = (min(hai_E_range[1], hc_E_range[1]), max(hai_E_range[2], hc_E_range[2]))

    return (
        hai_buffer=hai_buffer,
        hai_I=hai_I,
        hai_E=hai_E,
        hai_E_range=hai_E_range,
        hc_buffer=hc_buffer,
        hc_I=hc_I,
        hc_E=hc_E,
        hc_E_range=hc_E_range,
        diff_abs=diff_abs,
        diff_rel=diff_rel,
        lims_err=lims_err,
    )
end

function _print_comparison(stats)
    print(
        """
HAI:
  I ....... = $(stats.hai_I)
  E ....... = $(@sprintf("%.6e", stats.hai_E))
  #boxes .. = $(length(stats.hai_buffer))
  err_range = ($(@sprintf("%.2e", stats.hai_E_range[1])), $(@sprintf("%.2e", stats.hai_E_range[2])))

HC:
  I ....... = $(stats.hc_I)
  E ....... = $(@sprintf("%.6e", stats.hc_E))
  #boxes .. = $(length(stats.hc_buffer))
  err_range = ($(@sprintf("%.2e", stats.hc_E_range[1])), $(@sprintf("%.2e", stats.hc_E_range[2])))

diff_abs = $(@sprintf("%.2e", stats.diff_abs))
diff_rel = $(@sprintf("%.2e", stats.diff_rel))
        """,
    )
    return nothing
end

function draw_2d_boxes!(ax, lows, upps, cvalues, crange; cmap=:viridis, draw_edges=false)
    n = length(cvalues)
    rects = Vector{Rect2f}(undef, n)

    @inbounds for i in eachindex(cvalues)
        low = Point2f(lows[i]...)
        high = Point2f(upps[i]...)
        rects[i] = Rect2f(low, Vec2f(high .- low))
    end

    poly!(
        ax,
        rects;
        color=cvalues,
        colormap=cmap,
        colorrange=crange,
        colorscale=log10,
        transparency=true,
        strokecolor=(draw_edges ? :black : :transparent),
        strokewidth=(draw_edges ? 1 : 0),
    )

    return nothing
end

function plot_buffer_2d(f, rtol::Real, draw_edges::Bool=false)
    a, b = (-1, -1), (1, 1)
    stats = _compute_buffers_and_stats(f, a, b, rtol)
    _print_comparison(stats)

    fig = Figure()
    axs = (
        Axis(fig[1, 1]; title="HAdaptiveIntegration"), Axis(fig[1, 3]; title="HCubature")
    )

    for ax in axs
        xlims!(ax, -1, 1)
        ylims!(ax, -1, 1)
        ax.aspect = DataAspect()
    end

    hai_a = map(x -> x[1].corners[1], stats.hai_buffer.valtree)
    hai_b = map(x -> x[1].corners[2], stats.hai_buffer.valtree)
    hai_E = map(x -> x[3], stats.hai_buffer.valtree)
    draw_2d_boxes!(axs[1], hai_a, hai_b, hai_E, stats.lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 2]; colormap=:viridis, limits=stats.lims_err, scale=log10)

    hc_a = map(x -> x.a, stats.hc_buffer.valtree)
    hc_b = map(x -> x.b, stats.hc_buffer.valtree)
    hc_E = map(x -> x.E, stats.hc_buffer.valtree)
    draw_2d_boxes!(axs[2], hc_a, hc_b, hc_E, stats.lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 4]; colormap=:viridis, limits=stats.lims_err, scale=log10)

    return fig
end

function draw_3d_boxes!(
    ax, lows, upps, cvalues, crange; cmap=:viridis, alpha=0.2, draw_edges=false
)
    n = length(cvalues)
    centers = Vector{Point3f}(undef, n)
    sizes = Vector{Vec3f}(undef, n)

    @inbounds for i in eachindex(cvalues)
        low = Point3f(lows[i]...)
        high = Point3f(upps[i]...)
        centers[i] = Point3f((low .+ high) ./ 2)
        sizes[i] = Vec3f(high .- low)
    end

    marker = Rect3f(Point3f(-0.5f0, -0.5f0, -0.5f0), Vec3f(1.0f0, 1.0f0, 1.0f0))

    meshscatter!(
        ax,
        centers;
        marker=marker,
        markersize=sizes,
        color=cvalues,
        colormap=cmap,
        colorrange=crange,
        colorscale=log10,
        transparency=true,
        alpha=alpha,
        shading=NoShading,
    )

    if draw_edges
        @inbounds for i in eachindex(cvalues)
            wireframe!(
                ax, Rect3f(centers[i] - sizes[i] ./ 2, sizes[i]); color=:black, linewidth=1
            )
        end
    end

    return nothing
end

function plot_buffer_3d(f, rtol::Real, draw_edges::Bool=false)
    a, b = (-1, -1, -1), (1, 1, 1)
    stats = _compute_buffers_and_stats(f, a, b, rtol)
    _print_comparison(stats)

    fig = Figure()
    axs = (
        Axis3(fig[1, 1]; title="HAdaptiveIntegration"), Axis3(fig[1, 3]; title="HCubature")
    )

    for ax in axs
        xlims!(ax, -1, 1)
        ylims!(ax, -1, 1)
        zlims!(ax, -1, 1)
        ax.aspect = (1, 1, 1)
    end

    hai_a = map(x -> x[1].corners[1], stats.hai_buffer.valtree)
    hai_b = map(x -> x[1].corners[2], stats.hai_buffer.valtree)
    hai_E = map(x -> x[3], stats.hai_buffer.valtree)
    draw_3d_boxes!(axs[1], hai_a, hai_b, hai_E, stats.lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 2]; colormap=:viridis, limits=stats.lims_err, scale=log10)

    hc_a = map(x -> x.a, stats.hc_buffer.valtree)
    hc_b = map(x -> x.b, stats.hc_buffer.valtree)
    hc_E = map(x -> x.E, stats.hc_buffer.valtree)
    draw_3d_boxes!(axs[2], hc_a, hc_b, hc_E, stats.lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 4]; colormap=:viridis, limits=stats.lims_err, scale=log10)

    return fig
end
