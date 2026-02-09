"""
    struct Orthotope{D,T} <: AbstractDomain{D,T}

Axes-aligned Orthotope in `D` dimensions, with element type `T`, given by two points
`low_corner` and `high_corner`. Note that, we must have `low_corner .≤ high_corner`.

## Fields:
- `corners::SVector{2,SVector{D,T}}`: `corners[1]` is the low corner and `corners[2]` is the
   high corner.

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
        @assert length(high_corner) == D "`low_corner` and `high_corner` must have the \
        same length."
    else
        @assert length(low_corner) == D "low_corner must have length $D."
        @assert length(high_corner) == D "high_corner must have length $D."
    end
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner)) "must have `low_corner \
.≤ high_corner`."

    return Orthotope(SVector(SVector{D,T}(low_corner), SVector{D,T}(high_corner)))
end

function Orthotope(low_corner, high_corner, D::Union{Int,Nothing}=nothing)
    if isnothing(D)
        D = length(low_corner)
        @assert length(high_corner) == D "`low_corner` and `high_corner` must have the \
        same length."
    else
        @assert length(low_corner) == D "low_corner must have length $D."
        @assert length(high_corner) == D "high_corner must have length $D."
    end
    @assert all(a ≤ b for (a, b) in zip(low_corner, high_corner)) "must have `low_corner \
.≤ high_corner`."

    return Orthotope(SVector(float(SVector{D}(low_corner)), float(SVector{D}(high_corner))))
end

"""
    reference_orthotope(D::Int, T=float(Int))

Return the reference `D`-dimensional orthotope `[0, 1]ᴰ` with element type `T`.
"""
function reference_orthotope(::Val{D}, (::Type{T})=float(Int)) where {D,T}
    return Orthotope(SVector{2}(zeros(SVector{D,T}), ones(SVector{D,T})))
end

function map_from_reference(h::Orthotope{D,T}) where {D,T}
    diff = h.corners[2] - h.corners[1]
    return (u -> h.corners[1] .+ u .* diff, prod(diff))
end

function map_to_reference(h::Orthotope{D,T}) where {D,T}
    diff = h.corners[2] - h.corners[1]
    @assert all(x -> x > √eps(float(T)), diff) "degenerate $D-dimensional Orthotope: must \
have `high_corner .> low_corner`."

    return p -> (p - h.corners[1]) ./ diff
end

"""
    subdivide_orthotope(h::Orthotope)

Subdivide the `D`-orthotope `h` into 2ᴰ smaller orthotopes by splitting each dimension at
its midpoint.
"""
function subdivide_orthotope(h::Orthotope{D,T}) where {D,T}
    a = h.corners[1]
    b = h.corners[2]
    m = (a + b) / 2

    return ntuple(Val(2^D)) do k
        aₖ = MVector{D,T}(undef)
        bₖ = MVector{D,T}(undef)
        for j in 1:D
            # `k-1` encodes, in binary, which half (lower or upper) is selected for each
            # dimension `j`. Bit `j` of `k-1` indicates whether the sub-orthotope uses the
            # lower half along `j`.
            if isodd((k - 1) >> (j - 1))
                aₖ[j] = a[j]
                bₖ[j] = m[j]
            else
                aₖ[j] = m[j]
                bₖ[j] = b[j]
            end
        end
        return Orthotope(SVector(SVector{D}(aₖ), SVector{D}(bₖ)))
    end
end
