"""
    const Cuboid{T} = Orthotope{3,T}

Alias for a 3-dimensional [`Orthotope`](@ref) of value type `T`.

## Constructors:
- `Cuboid(low_corner, high_corner)`
- `Cuboid{T}(low_corner, high_corner)`
"""
const Cuboid{T} = Orthotope{3,T}

Cuboid{T}(low_corner, high_corner) where {T} = Orthotope{T}(low_corner, high_corner, 3)
Cuboid(low_corner, high_corner) = Orthotope{float(Int)}(low_corner, high_corner, 3)

"""
    subdivide_cuboid(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints
of the edges.
"""
function subdivide_cuboid(c::Cuboid{T}) where {T}
    a, b = c.corners
    m = (a + b) / 2
    return (
        Cuboid{T}(SVector(a, m)),
        Cuboid{T}(SVector(SVector(m[1], a[2], a[3]), SVector(b[1], m[2], m[3]))),
        Cuboid{T}(SVector(SVector(a[1], m[2], a[3]), SVector(m[1], b[2], m[3]))),
        Cuboid{T}(SVector(SVector(m[1], m[2], a[3]), SVector(b[1], b[2], m[3]))),
        Cuboid{T}(SVector(SVector(a[1], a[2], m[3]), SVector(m[1], m[2], b[3]))),
        Cuboid{T}(SVector(SVector(m[1], a[2], m[3]), SVector(b[1], m[2], b[3]))),
        Cuboid{T}(SVector(SVector(a[1], m[2], m[3]), SVector(m[1], b[2], b[3]))),
        Cuboid{T}(SVector(m, b)),
    )
end
