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

function integral_monomial_simplex(dim::Int, kmax::Int)
    base = Matrix(I, dim, dim)

    values = Dict(Tuple(0 for _ in 1:dim) => 1//factorial(dim))
    for k in 1:kmax
    end

    return values
end

function integral_monomial_orthotope() end
