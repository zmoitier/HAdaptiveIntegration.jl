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

function Orthotope{T}(low_corner, high_corner) where {T}
    @assert (length(low_corner) == length(high_corner)) "corners must have the same length."
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner)) "low corner must be less than or equal to high corner."

    D = length(low_corner)
    return Orthotope(SVector{D,T}(low_corner), SVector{D,T}(high_corner))
end
function Orthotope(low_corner, high_corner)
    return Orthotope{promote_to_float(low_corner, high_corner)}(low_corner, high_corner)
end

"""
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
