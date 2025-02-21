"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment{T}) where {T}
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment(a, m), Segment(m, b))
end

"""
    subdivide_triangle2(s::Triangle)

Divide the triangle `t` into two triangles by connecting the first point of `t` to the midpoints of
the two other points.
"""
function subdivide_triangle2(t::Triangle)
    a, b, c = t.points
    bc = (b + c) / 2
    return (Triangle(bc, a, b), Triangle(bc, c, a))
end

"""
    subdivide_triangle4(t::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle4(t::Triangle)
    a, b, c = t.points
    ab = (a + b) / 2
    ac = (c + a) / 2
    bc = (b + c) / 2
    return (
        Triangle(a, ab, ac), Triangle(b, bc, ab), Triangle(c, ac, bc), Triangle(ab, bc, ac)
    )
end

"""
    subdivide_rectangle4(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the midpoints
of the edges.
"""
function subdivide_rectangle4(r::Rectangle)
    lc = r.low_corner
    hc = r.high_corner
    m = (lc + hc) / 2
    return (
        Rectangle(r.low_corner, m),
        Rectangle(SVector{2}(m[1], lc[2]), SVector{2}(hc[1], m[2])),
        Rectangle(SVector{2}(lc[1], m[2]), SVector{2}(m[1], hc[2])),
        Rectangle(m, r.high_corner),
    )
end

"""
    subdivide_tetrahedron8(t::Tetrahedron)

Divide the tetrahedron `t` into eight tetrahedra by connecting the midpoints of the edges.
"""
function subdivide_tetrahedron8(t::Tetrahedron)
    a, b, c, d = t.points
    ab = (a + b) / 2
    ac = (a + c) / 2
    ad = (a + d) / 2
    bc = (b + c) / 2
    bd = (b + d) / 2
    cd = (c + d) / 2
    return [
        # (1/2)-tetrahedron on each vertices
        Tetrahedron(a, ab, ac, ad),
        Tetrahedron(ab, b, bc, bd),
        Tetrahedron(ac, bc, c, cd),
        Tetrahedron(ad, bd, cd, d),
        # octahedron splitting
        Tetrahedron(ab, bc, bd, ad),
        Tetrahedron(ac, bc, cd, ad),
        Tetrahedron(ad, cd, bd, bc),
        Tetrahedron(bc, ac, ab, ad),
    ]
end

"""
    subdivide_cuboid8(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints of the
edges.
"""
function subdivide_cuboid8(c::Cuboid)
    xl, yl, zl = c.low_corner
    xu, yu, zu = c.high_corner
    xm, ym, zm = (xl + xu) / 2, (yl + yu) / 2, (zl + zu) / 2
    return (
        Cuboid((xl, yl), (xm, ym), (zl, zm)),
        Cuboid((xm, yl), (xu, ym), (zl, zm)),
        Cuboid((xl, ym), (xm, yu), (zl, zm)),
        Cuboid((xm, ym), (xu, yu), (zl, zm)),
        Cuboid((xl, yl), (xm, ym), (zm, zu)),
        Cuboid((xm, yl), (xu, ym), (zm, zu)),
        Cuboid((xl, ym), (xm, yu), (zm, zu)),
        Cuboid((xm, ym), (xu, yu), (zm, zu)),
    )
end

default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cuboid) = subdivide_cuboid8
