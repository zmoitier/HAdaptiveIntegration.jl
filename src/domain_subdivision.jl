# TODO: move to utils
"""
    check_subdivision(
        domain::Domain{D,T},
        subdiv_algo;
        atol::T=zero(T),
        rtol::T=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {D,T<:Real}

Return `nothing` if the sum of the volume of the subdomain by the `subdiv_algo` is equal to the volume of the domain else throw an error.
"""
function check_subdivision(
    domain::Domain{D,T},
    subdiv_algo;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T<:Real}
    subdomains = subdiv_algo(domain)
    if isapprox(
        sum(abs_det_jacobian.(subdomains)), abs_det_jacobian(domain); atol=atol, rtol=rtol
    )
        @info "`$(Symbol(subdiv_algo))` pass volume test."
        return nothing
    else
        @error "`$(Symbol(subdiv_algo))` do not partition the domain."
    end
end

"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment{T}) where {T<:Real}
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment{T}(a, m), Segment{T}(m, b))
end

"""
    subdivide_segment3(s::Segment)

Divide the segment `s` into three segments of equal length.
"""
function subdivide_segment3(s::Segment{T}) where {T<:Real}
    a, b = s.low_corner, s.high_corner
    m1, m2 = (2 * a + b) / 3, (a + 2 * b) / 3
    return (Segment{T}(a, m1), Segment{T}(m1, m2), Segment{T}(m2, b))
end

"""
    subdivide_triangle2(s::Triangle)

Divide the triangle `t` into two triangles by connecting the first point of `t` to the midpoints of
the two other points.
"""
function subdivide_triangle2(t::Triangle{T}) where {T<:Real}
    a, b, c = t.points
    bc = (b + c) / 2
    return (Triangle{T}(bc, a, b), Triangle{T}(bc, c, a))
end

"""
    subdivide_triangle4(t::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle4(t::Triangle{T}) where {T<:Real}
    a, b, c = t.points
    ab = (a + b) / 2
    ac = (c + a) / 2
    bc = (b + c) / 2
    return (
        Triangle{T}(a, ab, ac),
        Triangle{T}(b, bc, ab),
        Triangle{T}(c, ac, bc),
        Triangle{T}(ab, bc, ac),
    )
end

"""
    subdivide_rectangle4(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the midpoints
of the edges.
"""
function subdivide_rectangle4(r::Rectangle{T}) where {T<:Real}
    a, b = r.low_corner, r.high_corner
    m = (a + b) / 2
    return (
        Rectangle{T}(a, m),
        Rectangle{T}(SVector(m[1], a[2]), SVector(b[1], m[2])),
        Rectangle{T}(SVector(a[1], m[2]), SVector(m[1], b[2])),
        Rectangle{T}(m, b),
    )
end

"""
    subdivide_tetrahedron8(t::Tetrahedron)

Divide the tetrahedron `t` into eight tetrahedra by connecting the midpoints of the edges.
"""
function subdivide_tetrahedron8(t::Tetrahedron{T}) where {T<:Real}
    a, b, c, d = t.points
    ab = (a + b) / 2
    ac = (a + c) / 2
    ad = (a + d) / 2
    bc = (b + c) / 2
    bd = (b + d) / 2
    cd = (c + d) / 2
    return [
        # (1/2)-tetrahedron on each vertices
        Tetrahedron{T}(a, ab, ac, ad),
        Tetrahedron{T}(ab, b, bc, bd),
        Tetrahedron{T}(ac, bc, c, cd),
        Tetrahedron{T}(ad, bd, cd, d),
        # octahedron splitting
        Tetrahedron{T}(ab, bc, bd, ad),
        Tetrahedron{T}(ac, bc, cd, ad),
        Tetrahedron{T}(ad, cd, bd, bc),
        Tetrahedron{T}(bc, ac, ab, ad),
    ]
end

"""
    subdivide_cuboid8(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints of the
edges.
"""
function subdivide_cuboid8(c::Cuboid{T}) where {T<:Real}
    a, b = c.low_corner, c.high_corner
    m = (a + b) / 2
    return (
        Cuboid{T}(a, m),
        Cuboid{T}(SVector(m[1], a[2], a[3]), SVector(b[1], m[2], m[3])),
        Cuboid{T}(SVector(a[1], m[2], a[3]), SVector(m[1], b[2], m[3])),
        Cuboid{T}(SVector(a[1], a[2], m[3]), SVector(m[1], m[2], b[3])),
        Cuboid{T}(SVector(a[1], m[2], m[3]), SVector(m[1], b[2], b[3])),
        Cuboid{T}(SVector(m[1], a[2], m[3]), SVector(b[1], m[2], b[3])),
        Cuboid{T}(SVector(m[1], m[2], a[3]), SVector(b[1], b[2], m[3])),
        Cuboid{T}(m, b),
    )
end

const LIST_SUBDIVISION_ALGO = [
    "segment" => ["subdivide_segment2", "subdivide_segment3"],
    "rectangle" => ["subdivide_rectangle4"],
    "triangle" => ["subdivide_triangle2", "subdivide_triangle4"],
    "cuboid" => ["subdivide_cuboid8"],
    "tetrahedron" => ["subdivide_tetrahedron8"],
]

function default_subdivision(d::Domain)
    @error "no default subdivision for $(typeof(d))."
end

default_subdivision(::Segment) = subdivide_segment2
default_subdivision(::Triangle) = subdivide_triangle4
default_subdivision(::Rectangle) = subdivide_rectangle4
default_subdivision(::Tetrahedron) = subdivide_tetrahedron8
default_subdivision(::Cuboid) = subdivide_cuboid8
