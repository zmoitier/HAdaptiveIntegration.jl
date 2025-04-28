"""
Alias for a 3-dimensional tetrahedron with 4 vertices of element type `T`.

## Constructors:
- `Tetrahedron(a, b, c, d)`
- `Tetrahedron{T}(a, b, c, d)`
"""
const Tetrahedron{T} = Simplex{3,T,4}

function Tetrahedron{T}(a, b, c, d) where {T}
    @assert length(a) == length(b) == length(c) == length(d) == 3 "all `vertices` must have length 2."
    return Simplex(
        SVector{4}(SVector{3,T}(a), SVector{3,T}(b), SVector{3,T}(c), SVector{3,T}(d))
    )
end
Tetrahedron(a, b, c, d) = Tetrahedron{promote_to_float(a, b, c, d)}(a, b, c, d)

"""
    subdivide_tetrahedron(t::Tetrahedron)

Divide the tetrahedron `t` into eight tetrahedra by connecting the midpoints of the edges.
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
        Tetrahedron{T}(SVector{4}(a, ab, ac, ad)),
        Tetrahedron{T}(SVector{4}(ab, b, bc, bd)),
        Tetrahedron{T}(SVector{4}(ac, bc, c, cd)),
        Tetrahedron{T}(SVector{4}(ad, bd, cd, d)),
        # octahedron splitting
        Tetrahedron{T}(SVector{4}(ab, ac, ad, bd)),
        Tetrahedron{T}(SVector{4}(ab, ac, bc, bd)),
        Tetrahedron{T}(SVector{4}(ac, ad, bd, cd)),
        Tetrahedron{T}(SVector{4}(ac, bc, bd, cd)),
    )
end
