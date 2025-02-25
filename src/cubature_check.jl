function check_order_simplex(
    nodes::SVector{N,SVector{D,T}}, weights::SVector{N,T}, order::Int
) where {N,D,T<:Real}
    return nothing
end

function check_order_orthotope(
    nodes::SVector{N,SVector{D,T}}, weights::SVector{N,T}, order::Int
) where {N,D,T<:Real}
    return nothing
end
