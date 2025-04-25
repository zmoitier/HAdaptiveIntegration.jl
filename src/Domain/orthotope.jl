"""
    struct Orthotope{D,T} <: AbstractDomain{D,T}

Axes-aligned Orthotope in `D` dimensions, with element type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `low_corner::SVector{D,T}`: the low corner.
- `high_corner::SVector{D,T}`: the high corner.

## Invariants (**not** check at construction):
- `low_corner .≤ high_corner`
"""
struct Orthotope{D,T} <: AbstractDomain{D,T}
    low_corner::SVector{D,T}
    high_corner::SVector{D,T}
end

"""
    orthotope(low_corner, high_corner)
    orthotope(T::DataType, low_corner, high_corner)

Return an axes-aligned orthotope in `D` dimensions, with element type `T`, given by two
points `low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function orthotope(T::DataType, low_corner, high_corner)
    @assert (length(low_corner) == length(high_corner))
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner))

    D = length(low_corner)
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end
function orthotope(low_corner, high_corner)
    return orthotope(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    reference_orthotope(D::Int)
    reference_orthotope(T::DataType, D::Int)

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with element type `T`.
"""
function reference_orthotope(T::DataType, D::Int)
    return Orthotope(zeros(SVector{D,T}), ones(SVector{D,T}))
end

function map_from_reference(h::Orthotope{D,T}) where {D,T}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

function abs_det_jac(h::Orthotope{D,T}) where {D,T}
    return prod(h.high_corner - h.low_corner)
end

function map_to_reference(h::Orthotope{D,T}) where {D,T}
    diff = h.high_corner - h.low_corner
    @assert all(x -> x > √eps(float(T)), diff) "degenerate $D-dimensional Orthotope: must have `high_corner .> low_corner`."

    return p -> (p - h.low_corner) ./ diff
end

"""
    subdivide_reference_orthotope(::Val{D}, ::Type{T}=Float64) where {D,T}

Like `subdivide_orthotope`, but operates on the reference orthotope. Since the output
depends only on the dimension `D`, and the type `T` used to represent coordinates, this
function is generated for each combination of `D` and `T`.
"""
@generated function subdivide_reference_orthotope(::Val{D}, (::Type{T})=Float64) where {D,T}
    a, b = zeros(SVector{D,T}), ones(SVector{D,T})
    m = SVector{D}(fill(T(1//2), D))

    sub_corners = Vector{NTuple{2,SVector{D,T}}}()
    for choices in Base.product([(true, false) for _ in 1:D]...)
        # Compute the low and high corners of the sub-orthotope
        low_corner = SVector{D,T}(cᵢ ? aᵢ : mᵢ for (cᵢ, aᵢ, mᵢ) in zip(choices, a, m))
        high_corner = SVector{D,T}(cᵢ ? mᵢ : bᵢ for (cᵢ, mᵢ, bᵢ) in zip(choices, m, b))
        push!(sub_corners, (low_corner, high_corner))
    end

    # Convert to an efficient format with known sizes
    static_orthotopes = ntuple(2^D) do i
        Orthotope(sub_corners[i][1], sub_corners[i][2])
    end

    return :($static_orthotopes)
end

"""
    subdivide_orthotope(h::Orthotope)

Subdivide the `D`-orthotope `h` into `2ᴰ` smaller orthotopes by splitting each dimension at
its midpoint.
"""
function subdivide_orthotope(h::Orthotope{D,T}) where {D,T}
    refs = subdivide_reference_orthotope(Val(D), T)
    f = map_from_reference(h)
    map(refs) do ref
        Orthotope(f(ref.low_corner), f(ref.high_corner))
    end
end

"""
    Segment{T} = Orthotope{1,T}

Alias for a 1-dimensional segment of element type `T`.
"""
const Segment{T} = Orthotope{1,T}

"""
    segment(xmin, xmax)
    segment(T::DataType, xmin, xmax)

Return a segment in 1 dimensions with element type `T` representing the interval
`[xmin, xmax]` with `xmin ≤ xmax`.
"""
function segment(T::DataType, xmin::R, xmax::S) where {R<:Real,S<:Real}
    return orthotope(T, (xmin,), (xmax,))
end
function segment(xmin::R, xmax::S) where {R<:Real,S<:Real}
    T = promote_type(typeof(xmin), typeof(xmax))
    return segment(float(T), xmin, xmax)
end

"""
    subdivide_segment2(s::Segment)

Divide the segment `s` into two segments of equal length.
"""
function subdivide_segment2(s::Segment{T}) where {T}
    a, b = s.low_corner, s.high_corner
    m = (a + b) / 2
    return (Segment{T}(a, m), Segment{T}(m, b))
end

"""
    Rectangle{T} = Orthotope{2,T}

Alias for a 2-dimensional rectangle of element type `T`.
"""
const Rectangle{T} = Orthotope{2,T}

"""
    rectangle(low_corner, high_corner)
    rectangle(T::DataType, low_corner, high_corner)

Return an axes-aligned rectangle, with element type `T`, given by two 2d-points `low_corner`
and `high_corner` with. Note that, we must have `low_corner .≤ high_corner`.
"""
function rectangle(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 2 "`low_corner` and `high_corner` must be 2d-vector."
    return orthotope(T, low_corner, high_corner)
end
function rectangle(low_corner, high_corner)
    return rectangle(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    subdivide_rectangle4(r::Rectangle)

Divide the rectangle `r` into four squares by connecting the center of the square to the
midpoints of the edges.
"""
function subdivide_rectangle4(r::Rectangle{T}) where {T}
    a, b = r.low_corner, r.high_corner
    m = (a + b) / 2
    return (
        Rectangle{T}(a, m),
        Rectangle{T}(SVector(m[1], a[2]), SVector(b[1], m[2])),
        Rectangle{T}(SVector(a[1], m[2]), SVector(m[1], b[2])),
        Rectangle{T}(m, b),
    )
end

"""
    Cuboid{T} = Orthotope{3,T}

Alias for a 3-dimensional cuboid of value type `T`.
"""
const Cuboid{T} = Orthotope{3,T}

"""
    cuboid(low_corner, high_corner)
    cuboid(T::DataType, low_corner, high_corner)

Return an axes-aligned cuboid, with element type `T`, given by two 3d-points `low_corner`
and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.
"""
function cuboid(T::DataType, low_corner, high_corner)
    @assert length(low_corner) == length(high_corner) == 3 "`low_corner` and `high_corner` must be 3d-vector."
    return orthotope(T, low_corner, high_corner)
end
function cuboid(low_corner, high_corner)
    return cuboid(promote_to_float(low_corner, high_corner), low_corner, high_corner)
end

"""
    subdivide_cuboid8(c::Cuboid)

Divide the cuboid `c` into 8 cuboid by connecting the center of the cuboid to the midpoints
of the edges.
"""
function subdivide_cuboid8(c::Cuboid{T}) where {T}
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
