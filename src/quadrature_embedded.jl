"""
    struct EmbeddedQuadrature{N,T}

An embedded quadrature rule consisting of a high order quadrature rule and a low order
quadrature rule.

## Fields:
- `nodes::Vector{SVector{N,T}}`: the quadrature nodes
- `weights_high::Vector{T}`: the quadrature weights for the high order
  quadrature
- `weights_low::Vector{T}`: the quadrature weights for the low order
  quadrature
- `label::String`: a label for the quadrature rule

Note that the low order quadrature uses `nodes[1:n]` as its nodes, where `n =
length(weights_low)`.
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
        @assert view(nodes_high, 1:length(nodes_low)) == nodes_low
        return new{N,T}(nodes_high, weights_high, weights_low, label)
    end
end

"""
    embedded_from_2quad(quad_low, quad_high, name, T = Float64)

Create an embedded quadrature rule from two quadrature rules. The nodes and
weights will be represented as `SVector{N,T}` and `T`, respectively.
"""
function embedded_from_2quad(quad_low, quad_high, name::String, T::DataType = Float64)
    _check_precision(quad_low, T)
    _check_precision(quad_high, T)
    N = length(quad_low.nodes[1])
    return EmbeddedQuadrature(
        SVector{N,T}.(quad_high.nodes),
        SVector{N,T}.(quad_low.nodes),
        T.(quad_high.weights),
        T.(quad_low.weights),
        name,
    )
end

function _check_precision(quad::Quadrature{<:Any,S}, T::DataType) where {S}
    if eps(T) < eps(S)
        error("requested precision $T cannot be achieved with quadrature rule of type {S}")
    end
end

function (quad::EmbeddedQuadrature)(fct::Function, domain, norm = LinearAlgebra.norm)
    mu = det_jac(domain)
    phi = map_from_ref(domain)
    x_ref = quad.nodes
    w_low = quad.weights_low
    w_high = quad.weights_high
    n_low, n_high = length(quad.weights_low), length(quad.weights_high)

    # assuming that nodes in quad_high are ordered so that the overlapping nodes
    # come first, add them up
    x = phi(x_ref[1])
    I_high = fct(x) * w_high[1]
    I_low = fct(x) * w_low[1]
    for i = 2:n_low
        x = phi(x_ref[i])
        v = fct(x)
        I_high += v * w_high[i]
        I_low += v * w_low[i]
    end

    # now compute the rest of the high order quadrature
    for i = n_low+1:n_high
        x = phi(x_ref[i])
        I_high += fct(x) * w_high[i]
    end

    # return the integral and the error estimate
    return mu * I_high, mu * norm(I_high - I_low)
end

# default quadrature rules

const PREDEFINED_QUADRATURES =
    Set(["segment-G7K15", "triangle-LaurieRadon", "square-CoolsHaegemans"])

"""
    EmbeddedQuadrature(; name::String)

Create an embedded quadrature rule with the given `name`.
"""
function EmbeddedQuadrature(; name::String, datatype::DataType = Float64)
    name in PREDEFINED_QUADRATURES || error(
        "Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_QUADRATURES, ", "))",
    )

    if name == "segment-G7K15"
        return embedded_from_2quad(
            SEGMENT_GAUSS_O13_N7,
            SEGMENT_KRONROD_O23_N15,
            "G7K15",
            datatype,
        )
    elseif name == "triangle-LaurieRadon"
        return embedded_from_2quad(
            TRIANGLE_RADON_O5_N7,
            TRIANGLE_LAURIE_O8_N19,
            "LaurieRadon",
            datatype,
        )
    elseif name == "square-CoolsHaegemans"
        return embedded_from_2quad(
            SQUARE_COOLS_HAEGEMANS_O7_N21,
            SQUARE_GAUSS_O9_N25,
            "CoolsHaegemans",
            datatype,
        )
    else
        error("Unknown rule.")
    end
end

function default_quadrature(::Segment{T}) where {T}
    return EmbeddedQuadrature(; name = "segment-G7K15", datatype = T)
end
function default_quadrature(::Triangle{T}) where {T}
    return EmbeddedQuadrature(; name = "triangle-LaurieRadon", datatype = T)
end
function default_quadrature(::Square{T}) where {T}
    return EmbeddedQuadrature(; name = "square-CoolsHaegemans", datatype = T)
end
