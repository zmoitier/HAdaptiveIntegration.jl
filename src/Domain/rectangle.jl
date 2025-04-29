"""
    const Rectangle{T} = Orthotope{2,T}

Alias for a 2-dimensional [`Orthotope`](@ref) of element type `T`.

## Constructors:
- `Rectangle(low_corner, high_corner)`
- `Rectangle{T}(low_corner, high_corner)`
"""
const Rectangle{T} = Orthotope{2,T}

function Rectangle{T}(low_corner, high_corner) where {T}
    _validate_invariant_orthotope(2, low_corner, high_corner)
    return Orthotope{T}(low_corner, high_corner)
end

function Rectangle(low_corner, high_corner)
    _validate_invariant_orthotope(2, low_corner, high_corner)
    return Orthotope(low_corner, high_corner)
end

"""
    subdivide_rectangle(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the
midpoints of the edges.
"""
function subdivide_rectangle(r::Rectangle{T}) where {T}
    a, b = r.corners
    m = (a + b) / 2
    return (
        Rectangle{T}(SVector(a, m)),
        Rectangle{T}(SVector(SVector(m[1], a[2]), SVector(b[1], m[2]))),
        Rectangle{T}(SVector(SVector(a[1], m[2]), SVector(m[1], b[2]))),
        Rectangle{T}(SVector(m, b)),
    )
end
