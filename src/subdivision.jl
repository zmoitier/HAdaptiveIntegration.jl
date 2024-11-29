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
    subdivide_triangle4(s::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the
edges.
"""
function subdivide_triangle4(t::Triangle)
    p1, p2, p3 = t.points
    p12 = (p1 + p2) / 2
    p23 = (p2 + p3) / 2
    p31 = (p3 + p1) / 2
    return (
        Triangle(p1, p12, p31),
        Triangle(p2, p23, p12),
        Triangle(p3, p31, p23),
        Triangle(p12, p23, p31),
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

default_subdivision(::Segment)  = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Square)   = subdivide_square4
