"""
Alias for a 3-dimensional cuboid of value type `T`.

## Constructors:
- `Cuboid(low_corner, high_corner)`
- `Cuboid{T}(low_corner, high_corner)`
"""
const Cuboid{T} = Orthotope{3,T}

function Cuboid{T}(low_corner, high_corner) where {T}
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner` must have length 3."
    return Orthotope(SVector{3,T}(low_corner), SVector{3,T}(high_corner))
end
function Cuboid(low_corner, high_corner)
    return Cuboid{promote_to_float(low_corner, high_corner)}(low_corner, high_corner)
end

"""
    subdivide_cuboid(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints
of the edges.
"""
function subdivide_cuboid(c::Cuboid{T}) where {T}
    a, b = c.low_corner, c.high_corner
    m = (a + b) / 2
    return (
        Cuboid{T}(a, m),
        Cuboid{T}(SVector(m[1], a[2], a[3]), SVector(b[1], m[2], m[3])),
        Cuboid{T}(SVector(a[1], m[2], a[3]), SVector(m[1], b[2], m[3])),
        Cuboid{T}(SVector(m[1], m[2], a[3]), SVector(b[1], b[2], m[3])),
        Cuboid{T}(SVector(a[1], a[2], m[3]), SVector(m[1], m[2], b[3])),
        Cuboid{T}(SVector(m[1], a[2], m[3]), SVector(b[1], m[2], b[3])),
        Cuboid{T}(SVector(a[1], m[2], m[3]), SVector(m[1], b[2], b[3])),
        Cuboid{T}(m, b),
    )
end
