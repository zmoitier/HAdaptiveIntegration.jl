"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment)
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment(a, m), Segment(m, b))
end

"""
    subdivide_triangle2(s::Triangle)

Divide the triangle `t` into two triangles by connecting the midpoints of the
edges.
+
|\
| \
| /\
|/  \
+----+
"""
function subdivide_triangle2(t::Triangle)
    a, b, c = t.points
    bc = (b + c) / 2
    return (Triangle(bc, a, b), Triangle(bc, c, a))
end

"""
    subdivide_triangle4(s::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the
edges.
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
    subdivide_square4(s::Square)

Divide the square `s` into four squares by connecting the center of the square
to the midpoints of the edges.
"""
function subdivide_square4(s::Square)
    xl, yl = s.low_corner
    xu, yu = s.high_corner
    xm, ym = (xl + xu) / 2, (yl + yu) / 2
    return (
        Square((xl, yl), (xm, ym)), # lower left
        Square((xm, yl), (xu, ym)), # lower right
        Square((xl, ym), (xm, yu)), # upper left
        Square((xm, ym), (xu, yu)), # upper right
    )
end

"""
    subdivide_tetrahedron8(s::Tetrahedron)

Divide the triangle `t` into four triangles by connecting the midpoints of the
edges.
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
    subdivide_cube8(c::Cube)

Divide the square `s` into four squares by connecting the center of the square
to the midpoints of the edges.
"""
function subdivide_cube8(c::Cube)
    xl, yl, zl = c.low_corner
    xu, yu, zu = c.high_corner
    xm, ym, zm = (xl + xu) / 2, (yl + yu) / 2, (zl + zu) / 2
    return (
        Cube((xl, yl), (xm, ym), (zl, zm)),
        Cube((xm, yl), (xu, ym), (zl, zm)),
        Cube((xl, ym), (xm, yu), (zl, zm)),
        Cube((xm, ym), (xu, yu), (zl, zm)),
        Cube((xl, yl), (xm, ym), (zm, zu)),
        Cube((xm, yl), (xu, ym), (zm, zu)),
        Cube((xl, ym), (xm, yu), (zm, zu)),
        Cube((xm, ym), (xu, yu), (zm, zu)),
    )
end

default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Square) = subdivide_square4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cube) = subdivide_cube8
