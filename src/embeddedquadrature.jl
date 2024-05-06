"""
    struct EmbeddedQuadrature{N,T}

An embedded quadrature rule for the `N`-dimensional reference simplex defined by
the N+1 vertices `(0,...,0),(0,...,0,1),(0,...,0,1,0),(1,0,...,0)`. The low
order quadrature uses the first `n` `nodes`, where `n = length(weights_low)`.
"""
struct EmbeddedQuadrature{N,T}
    nodes::Vector{SVector{N,T}}
    weights_high::Vector{T}
    weights_low::Vector{T}
    label::String
    function EmbeddedQuadrature(
        nodes_high::Vector{SVector{N,T}},
        nodes_low::Vector{SVector{N,T}},
        weights_high::Vector{T},
        weights_low::Vector{T},
        label::String,
    ) where {N,T}
        @assert nodes_high[1:length(nodes_low)] == nodes_low
        return new{N,T}(nodes_high, weights_high, weights_low, label)
    end
end

const PREDEFINED_QUADRATURES = ["LaurieRadon"]

function EmbeddedQuadrature(; name::String)
    name in PREDEFINED_QUADRATURES || error(
        "Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_QUADRATURES, ", "))",
    )
    if name == "LaurieRadon"
        node_high    = TRIANGLE_L8N19[1]
        nodes_low    = TRIANGLE_R5N7[1]
        weights_high = TRIANGLE_L8N19[2]
        weights_low  = TRIANGLE_R5N7[2]
        return EmbeddedQuadrature(node_high, nodes_low, weights_high, weights_low, name)
    else
        error()
    end
end

function _integrate_with_error(
    f,
    t::Simplex{N},
    quad::EmbeddedQuadrature{N,T},
    norm = LinearAlgebra.norm,
) where {N,T}
    mu     = 2 * measure(t)
    phi    = parametrization(t)
    xref   = quad.nodes
    w_high = quad.weights_high
    w_low  = quad.weights_low
    nhigh  = length(w_high)
    nlow   = length(w_low)
    # assuming that nodes in quad_high are ordered so that the overlapping nodes
    # come first, add them up
    S      = Base.promote_op((x, w) -> f(x) * w, SVector{N,T}, T)
    I_high = zero(S)
    I_low  = zero(S)
    for i in 1:nlow
        x = phi(xref[i])
        v = f(x)
        I_high += v * w_high[i]
        I_low += v * w_low[i]
    end
    # now compute the rest of the high order quadrature
    for i in nlow+1:nhigh
        x = phi(xref[i])
        v = f(x)
        I_high += v * w_high[i]
    end
    return mu * I_high, mu * norm(I_high - I_low)
end

# default quadrature rule for simplices

#TODO: use type information when creating defaults (e.g. single-precision)
default_quadrature(::Triangle) = EmbeddedQuadrature(; name = "LaurieRadon")
