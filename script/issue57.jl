using BenchmarkTools
using GLMakie
using HAdaptiveIntegration: HAdaptiveIntegration as HAI
using HCubature: HCubature as HC
using Printf: @printf

function get_range_E(hai_buffer, hc_buffer)
    println("Range of E:")
    a, b = extrema(E for (_, _, E) in hai_buffer.valtree)
    @printf("HAI: %1.2e %1.2e\n", a, b)

    c, d = extrema(box.E for box in hc_buffer.valtree)
    @printf("HC:  %1.2e %1.2e\n", c, d)

    return (min(a, c), max(b, d))
end

function draw_boxes!(
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

function main(rtol::Real, draw_edges::Bool=false)
    # f((x₁, x₂, x₃),) = (√(1 - x₁^2) + √(1 - x₂^2) + √(1 - x₃^2))
    f(x) = sum(sqrt.(1 .- x .^ 2))

    hai_buffer = HAI.allocate_buffer(f, HAI.Cuboid((0, 0, 0), (1, 1, 1)))
    hai_I, hai_E = HAI.integrate(
        f, HAI.Cuboid((0, 0, 0), (1, 1, 1)); buffer=hai_buffer, rtol=rtol, maxsubdiv=10^12
    )

    hc_buffer = HC.hcubature_buffer(f, (0, 0, 0), (1, 1, 1))
    hc_I, hc_E = HC.hcubature(f, (0, 0, 0), (1, 1, 1); buffer=hc_buffer, rtol=rtol)

    println("buffer length:")
    println("HAI: ", length(hai_buffer.valtree))
    println("HC:  ", length(hc_buffer.valtree))

    lims_err = get_range_E(hai_buffer, hc_buffer)

    fig = Figure()
    axs = (
        Axis3(fig[1, 1]; title="HAdaptiveIntegration"), Axis3(fig[1, 3]; title="HCubature")
    )

    for ax in axs
        xlims!(ax, 0, 1)
        ylims!(ax, 0, 1)
        zlims!(ax, 0, 1)
        ax.aspect = (1, 1, 1)
    end

    hai_a = map(x -> x[1].corners[1], hai_buffer.valtree)
    hai_b = map(x -> x[1].corners[2], hai_buffer.valtree)
    hai_E = map(x -> x[3], hai_buffer.valtree)
    draw_boxes!(axs[1], hai_a, hai_b, hai_E, lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 2]; colormap=:viridis, limits=lims_err, scale=log10)

    hc_a = map(x -> x.a, hc_buffer.valtree)
    hc_b = map(x -> x.b, hc_buffer.valtree)
    hc_E = map(x -> x.E, hc_buffer.valtree)
    draw_boxes!(axs[2], hc_a, hc_b, hc_E, lims_err; draw_edges=draw_edges)
    Colorbar(fig[1, 4]; colormap=:viridis, limits=lims_err, scale=log10)

    return fig
end
