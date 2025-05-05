"""
    const Triangle{T} = Simplex{2,T,3}

Alias for a 2-dimensional [`Simplex`](@ref) with 3 vertices of value type `T`.

## Constructors:
- `Triangle(a, b, c)`
- `Triangle{T}(a, b, c)`
"""
const Triangle{T} = Simplex{2,T,3}

Triangle{T}(a, b, c) where {T} = Simplex{T}(a, b, c)
Triangle(a, b, c) = Simplex{float(Int)}(a, b, c)

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
        Triangle{T}(SVector(a, ab, ac)),
        Triangle{T}(SVector(b, bc, ab)),
        Triangle{T}(SVector(c, ac, bc)),
        Triangle{T}(SVector(ab, bc, ac)),
    )
end
