using CairoMakie
using CairoMakie.Makie.GeometryBasics
using HAdaptiveIntegration.Domain: Cuboid, Orthotope, Rectangle, Simplex, Tetrahedron, Triangle, reference_domain, subdivide_cuboid, subdivide_rectangle, subdivide_tetrahedron, subdivide_triangle
using LinearAlgebra
using StaticArrays

# Mix a color toward white by fraction f (f=1 keeps the color, f=0 is white).
# Lightens/desaturates without transparency.
function dim(color, f)
    c = f * RGBf(color) + (1 - f) * RGBf(1, 1, 1)
    return RGBAf(c.r, c.g, c.b, 1)
end

const COLORS2 = [dim(Makie.cgrad(:tab10, 10, categorical = true)[i], 0.35) for i in 1:8]
const COLORS3 = [dim(Makie.cgrad(:tab10, 10, categorical = true)[i], 0.65) for i in 1:8]

const EDGE_COLOR = :black

const TET_FACES = Makie.GLTriangleFace[
    Makie.GLTriangleFace(1, 3, 2),
    Makie.GLTriangleFace(1, 2, 4),
    Makie.GLTriangleFace(1, 3, 4),
    Makie.GLTriangleFace(2, 3, 4),
]

function barycenter(s::Simplex{D, T, N}) where {D, T, N}
    return sum(s.vertices) / N
end

function barycenter(h::Orthotope{D, T}) where {D, T}
    return (h.corners[1] + h.corners[2]) / 2
end

function translate(s::Simplex{D, T, N}, t::SVector{D, T}) where {D, T, N}
    return Simplex{D, T, N}(map(v -> v + t, s.vertices))
end

function translate(h::Orthotope{D, T}, t::SVector{D, T}) where {D, T}
    return Orthotope{D, T}(map(c -> c + t, h.corners))
end

function explode_offset(sub, b, factor)
    return factor * (barycenter(sub) - b)
end

function view_direction(azimuth, elevation)
    return SVector{3, Float64}(cos(elevation) * cos(azimuth), cos(elevation) * sin(azimuth), sin(elevation))
end

function plot_subdivisions!(ax, domain, subdivide, draw!, offset, view_dir = nothing; colors = COLORS2)
    b = barycenter(domain)
    subs = collect(subdivide(domain))
    if view_dir !== nothing
        # Painter's algorithm: draw far-to-near so opaque meshes render in correct order.
        sort!(subs; by = s -> dot(barycenter(s) + explode_offset(s, b, offset), view_dir))
    end
    for (sub, color) in zip(subs, colors)
        draw!(ax, translate(sub, explode_offset(sub, b, offset)), color)
    end
    return nothing
end

function simplex_vertices(s::Simplex{D}) where {D}
    return [Point{D, Float32}(Tuple(v)) for v in s.vertices]
end

function rectangle_vertices(rect)
    a, b = rect.corners
    return [Point2f(a[1], a[2]), Point2f(b[1], a[2]), Point2f(b[1], b[2]), Point2f(a[1], b[2])]
end

function plot_triangle!(fig, col)
    ax = Axis(fig[2, col]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    plot_subdivisions!(ax, reference_domain(Triangle), subdivide_triangle, draw_triangle!, 0.1)
    return nothing
end

function draw_triangle!(ax, tri, color)
    poly!(ax, collect(simplex_vertices(tri)); color = color, strokecolor = EDGE_COLOR, strokewidth = 1)
    return nothing
end

function plot_rectangle!(fig, col)
    ax = Axis(fig[2, col]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    plot_subdivisions!(ax, reference_domain(Rectangle), subdivide_rectangle, draw_rectangle!, 0.1)
    return nothing
end

function draw_rectangle!(ax, rect, color)
    poly!(ax, collect(rectangle_vertices(rect)); color = color, strokecolor = EDGE_COLOR, strokewidth = 1)
    return nothing
end

function plot_tetrahedron!(fig, col)
    azimuth, elevation = π / 4, π / 8
    ax = Axis3(
        fig[2, col]; aspect = :data, xlabel = "", ylabel = "", zlabel = "",
        azimuth = azimuth, elevation = elevation
    )
    hidedecorations!(ax)
    hidespines!(ax)
    plot_subdivisions!(ax, reference_domain(Tetrahedron), subdivide_tetrahedron, draw_tetrahedron!, 0.5, view_direction(azimuth, elevation); colors = COLORS3)
    return nothing
end

function draw_tetrahedron!(ax, tet, color)
    msh = GeometryBasics.Mesh(collect(simplex_vertices(tet)), TET_FACES)
    mesh!(ax, msh; color = color)
    wireframe!(ax, msh; color = EDGE_COLOR)
    return nothing
end

function plot_cuboid!(fig, col)
    azimuth, elevation = 2 * π / 3, π / 4
    ax = Axis3(
        fig[2, col]; aspect = :data, xlabel = "", ylabel = "", zlabel = "",
        azimuth = azimuth, elevation = elevation
    )
    hidedecorations!(ax)
    hidespines!(ax)
    plot_subdivisions!(ax, reference_domain(Cuboid), subdivide_cuboid, draw_cuboid!, 0.5, view_direction(azimuth, elevation); colors = COLORS3)
    return nothing
end

function draw_cuboid!(ax, cub, color)
    a, b = cub.corners
    rect = Rect3f(Point3f(a[1], a[2], a[3]), Vec3f(b[1] - a[1], b[2] - a[2], b[3] - a[3]))
    mesh!(ax, rect; color = color)
    wireframe!(ax, rect; color = EDGE_COLOR)
    return nothing
end

function main()
    fig = Figure(size = (1200, 400))

    # Shared title row so 2D (Axis) and 3D (Axis3) titles share one baseline.
    titles = ("Triangle", "Rectangle", "Tetrahedron", "Cuboid")
    for (col, title) in enumerate(titles)
        Label(fig[1, col], title; fontsize = 24, font = :bold, halign = :center, tellheight = true)
        colsize!(fig.layout, col, Relative(1 / 4))
    end

    plot_triangle!(fig, 1)
    plot_rectangle!(fig, 2)
    plot_tetrahedron!(fig, 3)
    plot_cuboid!(fig, 4)

    rowgap!(fig.layout, 1, 0)

    save("./paper/image/subdivision.pdf", fig)

    return nothing
end

main()
