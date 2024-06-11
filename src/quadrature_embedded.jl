"""
    struct EmbeddedQuadrature{N,T}

An embedded quadrature rule for the reference N-simplex. The low
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

"""
    size(quad::EmbeddedQuadrature)

Return the number of nodes in the low and high order quadrature.
"""
size(quad::EmbeddedQuadrature) = (length(quad.weights_low), length(quad.weights_high))

"""
    embedded_from_2quad(quad_low, quad_high, name)

Create an embedded quadrature rule from two quadrature rules.
"""
function embedded_from_2quad(
    quad_low::Quadrature{N,T},
    quad_high::Quadrature{N,T},
    name::String,
) where {N,T}
    return EmbeddedQuadrature(
        quad_high.nodes,
        quad_low.nodes,
        quad_high.weights,
        quad_low.weights,
        name,
    )
end

"""
    (quad::EmbeddedQuadrature{N,T})(
    fct::Function,
    simplex::Simplex{N},
    norm = LinearAlgebra.norm,
) where {N,T}

Compute the integral of `fct` over the simplex using the embedded quadrature.
"""
function (quad::EmbeddedQuadrature{N,T})(
    fct::Function,
    simplex::Simplex{N},
    norm = LinearAlgebra.norm,
) where {N,T}
    mu            = det_jac(simplex)
    phi           = map_from_ref(simplex)
    x_ref         = quad.nodes
    w_low         = quad.weights_low
    w_high        = quad.weights_high
    n_low, n_high = size(quad)

    # assuming that nodes in quad_high are ordered so that the overlapping nodes
    # come first, add them up
    x      = phi(x_ref[1])
    I_high = fct(x) * w_high[1]
    I_low  = fct(x) * w_low[1]
    for i in 2:n_low
        x = phi(x_ref[i])
        v = fct(x)
        I_high += v * w_high[i]
        I_low += v * w_low[i]
    end

    # now compute the rest of the high order quadrature
    for i in n_low+1:n_high
        x = phi(x_ref[i])
        I_high += fct(x) * w_high[i]
    end

    return mu * I_high, mu * norm(I_high - I_low)
end

# default quadrature rule for simplices

const PREDEFINED_QUADRATURES = Set(["segment-G7K15", "triangle-LaurieRadon"])

"""
    EmbeddedQuadrature(; name::String)

Create an embedded quadrature rule with the given `name`.
"""
function EmbeddedQuadrature(; name::String)
    name in PREDEFINED_QUADRATURES || error(
        "Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_QUADRATURES, ", "))",
    )

    if name == "segment-G7K15"
        return embedded_from_2quad(SEGMENT_G7, SEGMENT_K15, "G7K15")

    elseif name == "triangle-LaurieRadon"
        return embedded_from_2quad(TRIANGLE_R5N7, TRIANGLE_L8N19, "LaurieRadon")

    else
        error("Unknown rule.")
    end
end

#TODO: use type information when creating defaults (e.g. single-precision)
default_quadrature(::Segment) = EmbeddedQuadrature(; name = "segment-G7K15")
default_quadrature(::Triangle) = EmbeddedQuadrature(; name = "triangle-LaurieRadon")
