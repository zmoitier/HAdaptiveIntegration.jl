"""
    const Triangle{T} = Simplex{2,T,3}

Alias for a 2-dimensional [`Simplex`](@ref) with 3 vertices of value type `T`.

## Constructors:
- `Triangle(a, b, c)`
- `Triangle{T}(a, b, c)`
"""
const Triangle{T} = Simplex{2,T,3}

Triangle{T}(a, b, c) where {T} = Simplex{T}(a, b, c)
Triangle(a, b, c) = Simplex(a, b, c)

"""
    subdivide_triangle(t::Triangle)

Divide the triangle `t` into 4 triangles by connecting the midpoints of the edges.
"""
function subdivide_triangle(t::Triangle{T}) where {T}
    a, b, c = t.vertices
    ab = (a + b) / 2
    ac = (a + c) / 2
    bc = (b + c) / 2
    return (
        # (1/2)-triangle on each vertices
        Triangle{T}(SVector(a, ab, ac)),
        Triangle{T}(SVector(ab, b, bc)),
        Triangle{T}(SVector(ac, bc, c)),
        # middle triangle
        Triangle{T}(SVector(bc, ac, ab)),
    )
end
