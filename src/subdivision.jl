function segment_subdivide2(s::Segment{T}) where {T}
    a, b = s.points
    m = (a + b) / 2
    return (Segment(a, m), Segment(m, b))
end

function triangle_subdivide4(t::Triangle{T}) where {T}
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

const PREDEFINED_SUBDIVIDE = Set(["segment-2", "triangle-4"])

function subdivide(; name::String)
    name in PREDEFINED_SUBDIVIDE || error(
        "Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_SUBDIVIDE, ", "))",
    )

    if name == "segment-2"
        return segment_subdivide2

    elseif name == "triangle-4"
        return triangle_subdivide4

    else
        error("Unknown subdivide.")
    end
end

default_subdivision(::Segment)  = subdivide(; name = "segment-2")
default_subdivision(::Triangle) = subdivide(; name = "triangle-4")
