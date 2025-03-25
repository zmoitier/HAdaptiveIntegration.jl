using Base.Iterators: partition
using StaticArrays

struct AffineMap{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
end

function affine_map(A::AbstractMatrix, b::AbstractVector)
    n, m = size(A)
    @assert m == length(b)

    T = promote_type(eltype(A), eltype(b))
    return AffineMap(SMatrix{n,m,T}(A), SVector{m,T}(b))
end

function affine_map(A::AbstractMatrix)
    n, m = size(A)
    return AffineMap(SMatrix{n,m}(A), zero(SVector{m,eltype(A)}))
end

function (Φ::AffineMap)(x)
    return Φ.A * x + Φ.b
end

function (Φ::AffineMap{1,T,1})(x) where {T}
    return Φ.A[1] * x + Φ.b[1]
end

function assemble(
    order_high::Int,
    order_low::Int,
    points::Tuple{Vector{AffineMap{D,T,N}},String,String,String}...,
) where {D,T,N}
    orbit_maps = Vector{Vector{AffineMap{D,T,N}}}()
    nodes = Vector{SVector{D,T}}()
    weights_high = Vector{T}()
    weights_low = Vector{T}()

    for (maps, point, wh, wl) in points
        push!(orbit_maps, maps)
        push!(nodes, SVector{1}(parse(T, point)))
        push!(weights_high, parse(T, wh))
        if wl != ""
            push!(weights_low, parse(T, wl))
        end
    end

    return (
        orbit_maps=orbit_maps,
        nodes=nodes,
        weights_high=weights_high,
        weights_low=weights_low,
        order_high=order_high,
        order_low=order_low,
    )
end

# [ nodes[1]..., ..., nodes[end]..., weights_high..., weights_low... ]
function pack(
    nodes::Vector{SVector{D,T}}, weights_high::Vector{T}, weights_low::Vector{T}
) where {D,T}
    U = Vector{T}()
    for node in nodes
        append!(U, node)
    end
    append!(U, weights_high)
    append!(U, weights_low)
    return U, D, length(weights_high), length(weights_low)
end

function unpack(U::Vector{T}, orbit_maps::Vector{Vector{AffineMap{D,T,N}}}) where {D,T,N}
    H = length(orbit_maps)

    nodes = [SVector{D}(t) for t in partition(U[1:(H * D)], D)]
    weights_high = U[(H * D + 1):(H * D + H)]
    weights_low = U[(H * D + H + 1):end]

    return nodes, weights_high, weights_low
end
