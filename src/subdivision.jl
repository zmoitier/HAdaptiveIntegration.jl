function subdivide(s::Segment{T}) where {T}
    a, b = s.points
    m = (a + b) / 2
    return (Segment(a, m), Segment(m, b))
end

function subdivide(t::Triangle{T}) where {T}
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
