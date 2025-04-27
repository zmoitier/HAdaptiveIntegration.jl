"""
Alias for a 2-dimensional rectangle of element type `T`.

## Constructors:
- `Rectangle(low_corner, high_corner)`
- `Rectangle{T}(low_corner, high_corner)`
"""
const Rectangle{T} = Orthotope{2,T}

function Rectangle{T}(low_corner, high_corner) where {T}
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner` must have length 2."
    return Orthotope(SVector{2,T}(low_corner), SVector{2,T}(high_corner))
end
function Rectangle(low_corner, high_corner)
    return Rectangle{promote_to_float(low_corner, high_corner)}(low_corner, high_corner)
end

"""
    subdivide_rectangle(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the
midpoints of the edges.
"""
function subdivide_rectangle(r::Rectangle{T}) where {T}
    a, b = r.low_corner, r.high_corner
    m = (a + b) / 2
    return (
        Rectangle{T}(a, m),
        Rectangle{T}(SVector(m[1], a[2]), SVector(b[1], m[2])),
        Rectangle{T}(SVector(a[1], m[2]), SVector(m[1], b[2])),
        Rectangle{T}(m, b),
    )
end
