using GLMakie
using HAdaptiveIntegration.Domain:
    Simplex,
    Tetrahedron,
    reference_domain,
    subdivide_tetrahedron
using StaticArrays
using Printf

include("utils.jl")

# Mix a color toward white by fraction f (f = 1 keeps the color, f = 0 is white).
function dim(color, f)
    c = f * RGBf(color) + (1 - f) * RGBf(1, 1, 1)
    return RGBAf(c.r, c.g, c.b, 1)
end

const COLORS = [dim(Makie.cgrad(:tab10, 10, categorical = true)[i], 0.5) for i in reverse(1:8)]

function barycenter(s::Simplex{D, T, N}) where {D, T, N}
    return sum(s.vertices) / N
end

function translate(s::Simplex{D, T, N}, t::SVector{D, T}) where {D, T, N}
    return Simplex{D, T, N}(map(v -> v + t, s.vertices))
end

# Faces of a tetrahedron: one per 3-subset of the four vertices, in boundary order.
function tetrahedron_faces(s::Tetrahedron)
    vts = s.vertices
    return [
        [Point3f(vts[i]) for i in idx]
            for idx in ((1, 3, 2), (1, 2, 4), (1, 4, 3), (2, 3, 4))
    ]
end

# Triangulate a convex face given in boundary order: a triangle is used as-is, a quad
# is split along the (1, 3) diagonal. Flat vertex-index triplets, as taken by mesh!.
function triangulate_face(n::Int)
    n == 3 && return [1, 2, 3]
    n == 4 && return [1, 2, 3, 1, 3, 4]
    return error("cannot triangulate a face with ", n, " vertices")
end

function draw_faces!(ax::Axis3, faces)
    for (pts, color) in faces
        mesh!(ax, pts, triangulate_face(length(pts)); color = color)
        lines!(ax, push!(copy(pts), pts[1]); color = :black, linewidth = 0.5)
    end
    return nothing
end

# Wrap the azimuth into (-180°, 180°] for display.
function azimuth_deg(azimuth::Real)
    d = rad2deg(mod2pi(azimuth))
    return d > 180 ? d - 360 : d
end

function camera_title(azimuth::Real, elevation::Real)
    return @sprintf(
        "azimuth = %.1f°, elevation = %.1f°", azimuth_deg(azimuth), rad2deg(elevation)
    )
end

function main()
    set_theme!(paper_theme(16))
    fig = Figure(size = (700, 700))

    ref = reference_domain(Tetrahedron)
    b = barycenter(ref)

    faces = [
        (pts, color)
            for (dom, color) in zip(subdivide_tetrahedron(ref), COLORS)
            for pts in tetrahedron_faces(translate(dom, 0.5 * (barycenter(dom) - b)))
    ]

    ax = Axis3(fig[1, 1]; aspect = :data, azimuth = deg2rad(-120), elevation = deg2rad(-40))
    draw_faces!(ax, faces)

    # The camera observables are updated by the drag-rotate interaction, so the title
    # tracks the angles as the axis is rotated.
    ax.title[] = camera_title(ax.azimuth[], ax.elevation[])
    onany(ax.azimuth, ax.elevation) do azimuth, elevation
        ax.title[] = camera_title(azimuth, elevation)
    end

    return fig
end

main()
