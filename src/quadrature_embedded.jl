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

function _integrate_with_error(
    f,
    t::Simplex{N},
    quad::EmbeddedQuadrature{N,T},
    norm = LinearAlgebra.norm,
) where {N,T}
    mu     = 2 * measure(t)
    phi    = map_from_ref(t)
    xref   = quad.nodes
    w_high = quad.weights_high
    w_low  = quad.weights_low
    nhigh  = length(w_high)
    nlow   = length(w_low)
    # assuming that nodes in quad_high are ordered so that the overlapping nodes
    # come first, add them up
    x1     = phi(xref[1])
    I_high = f(x1) * w_high[1]
    I_low  = f(x1) * w_low[1]
    for i in 2:nlow
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

const PREDEFINED_QUADRATURES = ["segment-G7K15", "triangle-LaurieRadon"]

function EmbeddedQuadrature(; name::String)
    name in PREDEFINED_QUADRATURES || error(
        "Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_QUADRATURES, ", "))",
    )

    if name == "segment-G7K15"
        nodes_low    = SEGMENT_G7.nodes
        weights_low  = SEGMENT_G7.weights
        node_high    = SEGMENT_K15.nodes
        weights_high = SEGMENT_K15.weights
        return EmbeddedQuadrature(node_high, nodes_low, weights_high, weights_low, name)

    elseif name == "triangle-LaurieRadon"
        nodes_low    = TRIANGLE_R5N7.nodes
        weights_low  = TRIANGLE_R5N7.weights
        node_high    = TRIANGLE_L8N19.nodes
        weights_high = TRIANGLE_L8N19.weights
        return EmbeddedQuadrature(node_high, nodes_low, weights_high, weights_low, name)

    else
        error("Unknown rule.")
    end
end

#TODO: use type information when creating defaults (e.g. single-precision)
default_quadrature(::Segment) = EmbeddedQuadrature(; name = "segment-G7K15")
default_quadrature(::Triangle) = EmbeddedQuadrature(; name = "triangle-LaurieRadon")
