"""
    struct HyperRectangle{N,T} <: Domain{N,T}

Axis-aligned hyperrectangle in `N` dimensions given by
`low_corner::SVector{N,T}` and `high_corner::SVector{N,T}`.
"""
struct HyperRectangle{N,T} <: Domain{N,T}
    low_corner::SVector{N,T}
    high_corner::SVector{N,T}
end

const Segment{T} = HyperRectangle{1,T}
const Square{T} = HyperRectangle{2,T}
const Cube{T} = HyperRectangle{3,T}

# default types
Segment(args...) = Segment{Float64}(args...)
Square(args...) = Square{Float64}(args...)
Cube(args...) = Cube{Float64}(args...)

"""
    reference_hyperrectangle(N::Int, T::DataType)

Return the reference `HyperRectangle` in `N` dimensions, given by `[0, 1]^N`.
"""
function reference_hyperrectangle(N::Int, T::DataType)
    @assert N ≥ 1 "N = $N must be greater than 1."
    lc = @SVector zeros(T, N)
    hc = @SVector ones(T, N)
    return HyperRectange(lc, hc)
end

function map_from_ref(h::HyperRectangle{N,T})::Function where {N,T}
    return u -> h.low_corner + u .* (h.high_corner - h.low_corner)
end

"""
    det_jac(s::HyperRectangle)

The determinant of the Jacobian of the map from the reference simplex to the
physical simplex.
"""
function det_jac(s::HyperRectangle{T}) where {T}
    return prod(s.high_corner - s.low_corner)
end
