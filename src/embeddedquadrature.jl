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
end

const PREDEFINED_QUADRATURES = ["LaurieRadon"]

function EmbeddedQuadrature(; name::String)
    name in PREDEFINED_QUADRATURES || error("Unknown quadrature rule: $name. Options are: $(join(PREDEFINED_QUADRATURES, ", "))")
    if name == "LaurieRadon"
        nodes        = TRIANGLE_L8N19[1]
        weights_high = TRIANGLE_L8N19[2]
        weights_low  = TRIANGLE_R5N7[2]
        return EmbeddedQuadrature(nodes, weights_high, weights_low, name)
    else
        error()
    end
end

function integrate_with_error(
    f,
    t::Simplex{N},
    q::EmbeddedQuadrature{N},
    norm = LinearAlgebra.norm,
) where {N}
    μ      = measure(t)
    x̂     = q.nodes
    w_high = q.weights_high
    w_low  = q.weights_low
    I_high = f(first(x2)) * first(w1)
    I_low  = f(first(x2)) * first(w2)
    nhigh  = length(w_high)
    nlow   = length(w_low)
    # assuming that nodes in quad_high are ordered so that the overlapping nodes
    # come first, add them up
    for i in 2:nlow
        x = t(x̂[i])
        v = f(x)
        I_high += v * w_high[i]
        I_low += v * w_low[i]
    end
    # now compute the rest of the high order quadrature
    for i in nlow+1:nhigh
        x = t(x̂[i])
        v = f(x)
        I_high += v * w_high[i]
    end
    return μ * I_high, μ * norm(I2 - I1)
end
