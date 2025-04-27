"""
Alias for a 2-dimensional triangle with 3 vertices of value type `T`.

## Constructors:
- `Triangle(a, b, c)`
- `Triangle{T}(a, b, c)`
"""
const Triangle{T} = Simplex{2,T,3}

function Triangle{T}(a, b, c) where {T}
    @assert length(a) == length(b) == length(c) == 2 "all `vertices` must have length 2."
    return Simplex(SVector{3}(SVector{2,T}(a), SVector{2,T}(b), SVector{2,T}(c)))
end
Triangle(a, b, c) = Triangle{promote_to_float(a, b, c)}(a, b, c)

"""
    subdivide_triangle(t::Triangle)

Divide the triangle `t` into four triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle(t::Triangle{T}) where {T}
    a, b, c = t.vertices
    ab = (a + b) / 2
    ac = (c + a) / 2
    bc = (b + c) / 2
    return (
        Triangle{T}(SVector{3}(a, ab, ac)),
        Triangle{T}(SVector{3}(b, bc, ab)),
        Triangle{T}(SVector{3}(c, ac, bc)),
        Triangle{T}(SVector{3}(ab, bc, ac)),
    )
end
