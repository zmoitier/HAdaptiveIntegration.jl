"""
    struct Orthotope{D,T} <: AbstractDomain{D,T}

Axes-aligned Orthotope in `D` dimensions, with element type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `corners::SVector{2,SVector{D,T}}`: `corners[1]` = the low corner and `corners[2]` = the
   ]high corner.

## Invariants (**not** check at construction):
- `corners[1] .≤ corners[2]`

## Constructors:
- `Orthotope(low_corner, high_corner)`
- `Orthotope{T}(low_corner, high_corner)`
- `Orthotope(corners::SVector{2,SVector{D,T}})`
"""
struct Orthotope{D,T} <: AbstractDomain{D,T}
    corners::SVector{2,SVector{D,T}}
end

function Orthotope{T}(low_corner, high_corner, D::Union{Int,Nothing}=nothing) where {T}
    if isnothing(D)
        D = length(low_corner)
        @assert length(high_corner) == D "`low_corner` and `high_corner` must have the same length."
    else
        @assert length(low_corner) == D "low_corner must have length $D."
        @assert length(high_corner) == D "high_corner must have length $D."
    end
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner)) "must have `low_corner .≤ high_corner`."

    return Orthotope(SVector(SVector{D,T}(low_corner), SVector{D,T}(high_corner)))
end

function Orthotope(low_corner, high_corner, D::Union{Int,Nothing}=nothing)
    if isnothing(D)
        D = length(low_corner)
        @assert length(high_corner) == D "`low_corner` and `high_corner` must have the same length."
    else
        @assert length(low_corner) == D "low_corner must have length $D."
        @assert length(high_corner) == D "high_corner must have length $D."
    end
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner)) "must have `low_corner .≤ high_corner`."

    return Orthotope(SVector(float(SVector{D}(low_corner)), float(SVector{D}(high_corner))))
end

"""
    reference_orthotope(D::Int, T=float(Int))

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with element type `T`.
"""
function reference_orthotope(D::Int, (::Type{T})=float(Int)) where {T}
    return Orthotope(SVector{2}(zeros(SVector{D,T}), ones(SVector{D,T})))
end

function map_from_reference(h::Orthotope{D,T}) where {D,T}
    return u -> h.corners[1] + u .* (h.corners[2] - h.corners[1])
end

function abs_det_jac(h::Orthotope{D,T}) where {D,T}
    return prod(h.corners[2] - h.corners[1])
end

function map_to_reference(h::Orthotope{D,T}) where {D,T}
    diff = h.corners[2] - h.corners[1]
    @assert all(x -> x > √eps(float(T)), diff) "degenerate $D-dimensional Orthotope: must have `high_corner .> low_corner`."

    return p -> (p - h.corners[1]) ./ diff
end

"""
    subdivide_reference_orthotope(::Val{D}, ::Type{T}=float(Int))

Like `subdivide_orthotope`, but operates on the reference orthotope. Since the output
depends only on the dimension `D`, and the type `T` used to represent coordinates, this
function is generated for each combination of `D` and `T`.
"""
@generated function subdivide_reference_orthotope(
    ::Val{D}, (::Type{T})=float(Int)
) where {D,T}
    a, b = zeros(SVector{D,T}), ones(SVector{D,T})
    m = SVector{D}(fill(T(1//2), D))

    sub_corners = Vector{SVector{2,SVector{D,T}}}()
    for choices in Base.product([(true, false) for _ in 1:D]...)
        # Compute the low and high corners of the sub-orthotope
        low_corner = SVector{D,T}(cᵢ ? aᵢ : mᵢ for (cᵢ, aᵢ, mᵢ) in zip(choices, a, m))
        high_corner = SVector{D,T}(cᵢ ? mᵢ : bᵢ for (cᵢ, mᵢ, bᵢ) in zip(choices, m, b))
        push!(sub_corners, SVector(low_corner, high_corner))
    end

    # Convert to an efficient format with known sizes
    static_orthotopes = ntuple(2^D) do i
        Orthotope(sub_corners[i])
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
        Orthotope(f.(ref.corners))
    end
end
