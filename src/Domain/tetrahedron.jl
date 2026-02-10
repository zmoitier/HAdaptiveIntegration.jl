"""
    const Tetrahedron{T} = Simplex{3,T,4}

Alias for a 3-dimensional [`Simplex`](@ref) with 4 vertices of element type `T`.

## Constructors:
- `Tetrahedron(a, b, c, d)`
- `Tetrahedron{T}(a, b, c, d)`
"""
const Tetrahedron{T} = Simplex{3,T,4}

Tetrahedron{T}(a, b, c, d) where {T} = Simplex{T}(a, b, c, d)
Tetrahedron(a, b, c, d) = Simplex(a, b, c, d)

"""
    subdivide_tetrahedron(t::Tetrahedron)

Divide the tetrahedron `t` into 8 tetrahedra by connecting the midpoints of the edges.
"""
function subdivide_tetrahedron(t::Tetrahedron{T}) where {T}
    a, b, c, d = t.vertices
    ab = (a + b) / 2
    ac = (a + c) / 2
    ad = (a + d) / 2
    bc = (b + c) / 2
    bd = (b + d) / 2
    cd = (c + d) / 2
    return (
        # (1/2)-tetrahedron on each vertices
        Tetrahedron{T}(SVector(a, ab, ac, ad)),
        Tetrahedron{T}(SVector(ab, b, bc, bd)),
        Tetrahedron{T}(SVector(ac, bc, c, cd)),
        Tetrahedron{T}(SVector(ad, bd, cd, d)),
        # middle octahedron splitting
        Tetrahedron{T}(SVector(ab, ac, ad, bd)),
        Tetrahedron{T}(SVector(ab, ac, bc, bd)),
        Tetrahedron{T}(SVector(ac, ad, bd, cd)),
        Tetrahedron{T}(SVector(ac, bc, bd, cd)),
    )
end
