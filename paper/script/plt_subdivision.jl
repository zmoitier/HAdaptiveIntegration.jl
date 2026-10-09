using GLMakie
using HAdaptiveIntegration.Domain:
    Cuboid,
    Orthotope,
    Rectangle,
    Simplex,
    Tetrahedron,
    Triangle,
    reference_domain,
    subdivide_cuboid,
    subdivide_rectangle,
    subdivide_tetrahedron,
    subdivide_triangle
using StaticArrays

include("utils.jl")

# Mix a color toward white by fraction f (f = 1 keeps the color, f = 0 is white).
function dim(color, f)
    c = f * RGBf(color) + (1 - f) * RGBf(1, 1, 1)
    return RGBAf(c.r, c.g, c.b, 1)
end

const TAB10 = Makie.cgrad(:tab10, 10, categorical = true)
const COLORS2 = [dim(TAB10[i], 0.5) for i in 1:4]
const COLORS3 = [dim(TAB10[i], 0.5) for i in reverse(1:8)]

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

# Plain Vector (not SVector) so that Makie can stroke the polygon outline.
function simplex_vertices(s::Simplex{D}) where {D}
    return collect(Point{D, Float32}.(s.vertices))
end

function rectangle_vertices(rect)
    a, b = rect.corners
    return [
        Point2f(a[1], a[2]),
        Point2f(b[1], a[2]),
        Point2f(b[1], b[2]),
        Point2f(a[1], b[2]),
    ]
end

# Faces of a tetrahedron: one per 3-subset of the four vertices, in boundary order.
function tetrahedron_faces(s::Tetrahedron)
    vts = s.vertices
    return [
        [Point3f(vts[i]) for i in idx]
            for idx in ((1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4))
    ]
end

# Faces of a cuboid; the six axis-aligned boundary quads.
function cuboid_faces(c::Cuboid)
    lo, hi = c.corners
    vts = [Point3f(p) for p in Iterators.product(zip(lo, hi)...)]
    idxs = (
        [1, 2, 4, 3],  # x = lo[1]
        [5, 6, 8, 7],  # x = hi[1]
        [1, 2, 6, 5],  # y = lo[2]
        [3, 4, 8, 7],  # y = hi[2]
        [1, 3, 7, 5],  # z = lo[3]
        [2, 4, 8, 6],  # z = hi[3]
    )
    return [[vts[i] for i in idx] for idx in idxs]
end

# Triangulate a convex face given in boundary order: a triangle is used as-is, a quad
# is split along the (1, 3) diagonal. Flat vertex-index triplets, as taken by mesh!.
function triangulate_face(n::Int)
    n == 3 && return [1, 2, 3]
    n == 4 && return [1, 2, 3, 1, 3, 4]
    return error("cannot triangulate a face with ", n, " vertices")
end

function draw_faces!(ax::Axis, faces)
    for (pts, color) in faces
        poly!(ax, pts; color = color, strokecolor = :black, strokewidth = 0.5)
    end
    return nothing
end

function draw_faces!(ax::Axis3, faces)
    for (pts, color) in faces
        mesh!(ax, pts, triangulate_face(length(pts)); color = color)
        lines!(ax, push!(copy(pts), pts[1]); color = :black, linewidth = 0.5)
    end
    return nothing
end

# Subdivide the reference domain and push each subdomain away from its shared
# barycenter (exploded view), collecting (points, color) faces for drawing.
function exploded_faces(domain, subdivide, colors, offset, parts_of)
    ref = reference_domain(domain)
    b = barycenter(ref)
    return [
        (pts, color)
            for (dom, color) in zip(subdivide(ref), colors)
            for pts in parts_of(translate(dom, offset * (barycenter(dom) - b)))
    ]
end

function plot_2d!(fig::Figure, col::Int; domain, subdivide, vertices, offset)
    faces = exploded_faces(domain, subdivide, COLORS2, offset, dom -> (vertices(dom),))

    ax = Axis(fig[2, col]; aspect = DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    draw_faces!(ax, faces)
    return nothing
end

function plot_3d!(
        fig::Figure, col::Int;
        domain, subdivide, faces_of, offset, azimuth, elevation
    )
    faces = exploded_faces(domain, subdivide, COLORS3, offset, faces_of)

    ax = Axis3(
        fig[2, col]; aspect = :data,
        azimuth = azimuth, elevation = elevation, protrusions = 0
    )
    hidedecorations!(ax)
    hidespines!(ax)
    draw_faces!(ax, faces)
    return nothing
end

function main()
    fig_w, fig_h, font_size = figure_sizes(1, 1 / 3)
    set_theme!(paper_theme(font_size))
    fig = Figure(size = (fig_w, fig_h))

    # Shared title row so 2D (Axis) and 3D (Axis3) titles share one baseline.
    titles = ("Triangle", "Rectangle", "Tetrahedron", "Cuboid")
    for (col, title) in enumerate(titles)
        Label(
            fig[1, col],
            title;
            fontsize = font_size,
            font = :bold,
            halign = :center,
            tellheight = true,
        )
        colsize!(fig.layout, col, Relative(1 / 4))
    end

    plot_2d!(
        fig, 1; domain = Triangle, subdivide = subdivide_triangle,
        vertices = simplex_vertices, offset = 0.1
    )
    plot_2d!(
        fig, 2; domain = Rectangle, subdivide = subdivide_rectangle,
        vertices = rectangle_vertices, offset = 0.1
    )
    plot_3d!(
        fig, 3; domain = Tetrahedron, subdivide = subdivide_tetrahedron,
        faces_of = tetrahedron_faces, offset = 0.5,
        azimuth = deg2rad(-120), elevation = deg2rad(-40)
    )
    plot_3d!(
        fig, 4; domain = Cuboid, subdivide = subdivide_cuboid,
        faces_of = cuboid_faces, offset = 0.85,
        azimuth = deg2rad(-58), elevation = deg2rad(25)
    )

    rowgap!(fig.layout, 1, 0)

    save("./paper/image/subdivision.png", fig; px_per_unit = 8)
    return nothing
end

main()
