"""
    struct Quadrature{N,T}

A quadrature rule (nodes and weights) on the reference N-simplex.
"""
struct Quadrature{N,T}
    nodes::Vector{SVector{N,T}}
    weights::Vector{T}

    function Quadrature(nodes::Vector{SVector{N,T}}, weights::Vector{T}) where {N,T}
        @assert length(nodes) == length(weights) "Should have the same number of nodes as weights."

        int_1 = sum(weights)
        if factorial(N) * int_1 ≉ T(1)
            @warn "The constant is not exactly integrated. \
            Maybe check if the quadrature is defined on the \
            reference N-simplex." abs(int_1 * factorial(N) - 1)
        end

        return new{N,T}(nodes, weights)
    end
end

length(quad::Quadrature) = length(quad.weights)

"""
    (quad::Quadrature)(fct::Function)

Compute `sum(quad.weights .* fct.(quad.nodes))`.
"""
function (quad::Quadrature)(fct::Function)
    return sum(quad.weights .* fct.(quad.nodes))
end

"""
    map_quad(quad::Quadrature{N,T}, simplex::Simplex{N,T,Np1}) where {N,T,Np1}

Maps a quadrature on the simplex `simplex` to the quadrature on the reference N-simplex.
"""
function map_quad(quad::Quadrature{N,T}, simplex::Simplex{N,T,Np1}) where {N,T,Np1}
    Φ = map_to_ref(simplex)
    c = 1 / (det_jac(simplex) * factorial(N))

    return Quadrature(Φ.(quad.nodes), c .* quad.weights)
end
