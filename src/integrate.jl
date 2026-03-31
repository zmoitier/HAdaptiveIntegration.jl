"""
    integrate(
        fct,
        domain::AbstractDomain{D};
        rule::EmbeddedCubature{D}=default_rule(domain),
        subdiv_algo=default_subdivision(domain),
        norm=LinearAlgebra.norm,
        buffer=nothing,
        atol=nothing,
        rtol=nothing,
        maxsubdiv=2^(13 + D),
        callback=(I, E, nb_subdiv, buffer) -> nothing,
    ) where {D}

Adaptively integrate `fct` over `domain`.

Return `(I, E)` where `I` is the integral estimate and `E` is an error estimate computed
from an embedded cubature pair.

## Arguments
- `fct`: a function that maps `SVector{D,T}` to a value in `K`. The return type `K` must
  support addition and multiplication by scalars of type `typeof(one(T))`.
- `domain::AbstractDomain{D,T}`: the integration domain. Currently, we support
  [`Segment`](@ref), [`Triangle`](@ref), [`Rectangle`](@ref), [`Tetrahedron`](@ref),
  [`Cuboid`](@ref), `D`-dimensional [`Simplex`](@ref), and `D`-dimensional
  [`Orthotope`](@ref).

## Optional arguments
- `rule::EmbeddedCubature{D,T}=default_rule(domain)`: the embedded cubature rule. Each
  supported domain has a [`default_rule`](@ref).
- `subdiv_algo=default_subdivision(domain)`: the subdivision algorithm, each domain has a
  [`default_subdivision`](@ref).
- `norm=LinearAlgebra.norm`: norm used to estimate the error.
- `buffer=nothing`: optional heap used by the adaptive algorithm. Reusing a buffer from
  [`allocate_buffer`](@ref) can reduce allocations when calling `integrate` repeatedly.
- `atol=nothing`: absolute tolerance, if `nothing` is passed, it will be set to `zero(E)`.
- `rtol=nothing`: relative tolerance, if `nothing` is passed, it will be set to
  `sqrt(eps(one(E)))`.
- `maxsubdiv=2^(13 + D)`: maximum number of subdivisions.
- `callback=(I, E, nb_subdiv, buffer) -> nothing`: a callback function called for each
  estimate of `I` and `E`, including the initial estimate (`nb_subdiv=0`) and after each
  subdivision. The callback receives the current integral `I`, error estimate `E`, number
  of subdivisions `nb_subdiv`, and `buffer`.

## Notes
- Iteration stops when `E ≤ atol` or `E ≤ rtol * norm(I)` or `nb_subdiv == maxsubdiv`.
"""
function integrate(
    fct,
    domain::AbstractDomain{D};
    rule::EmbeddedCubature{D}=default_rule(domain),
    subdiv_algo=default_subdivision(domain),
    norm=LinearAlgebra.norm,
    buffer=nothing,
    atol=nothing,
    rtol=nothing,
    maxsubdiv=2^(13 + D),
    callback=(_, _, _, _) -> nothing,
) where {D}
    return _integrate(
        fct, domain, rule, subdiv_algo, buffer, norm, atol, rtol, maxsubdiv, callback
    )
end

@noinline function _integrate(
    fct::FCT,
    domain::DOM,
    rule::EmbeddedCubature{D},
    subdiv_algo,
    buffer,
    norm,
    atol,
    rtol,
    maxsubdiv,
    callback,
) where {D,FCT,DOM<:AbstractDomain{D}}
    I, E = rule(fct, domain, norm)

    # initialize or reset the buffer
    buffer = if isnothing(buffer)
        BinaryHeap{Tuple{DOM,typeof(I),typeof(E)}}(Base.Order.By(last, Base.Order.Reverse))
    else
        empty!(buffer)
    end
    push!(buffer, (domain, I, E))

    # set default tolerances if not provided
    εₐ = something(atol, zero(E))
    εᵣ = something(rtol, sqrt(eps(one(E))))

    nb_subdiv = 0
    while true
        callback(I, E, nb_subdiv, buffer)

        # check termination conditions
        if (E ≤ εₐ) || (E ≤ εᵣ * norm(I)) || (nb_subdiv == maxsubdiv)
            break
        end

        # subdivide the domain with the largest error
        domain, I_dom, E_dom = pop!(buffer)
        I -= I_dom
        E -= E_dom
        for child in subdiv_algo(domain)
            I_child, E_child = rule(fct, child, norm)
            I += I_child
            E += E_child
            push!(buffer, (child, I_child, E_child))
        end
        nb_subdiv += 1
    end

    if nb_subdiv == maxsubdiv
        @warn "maximum number of subdivisions reached `maxsubdiv=$maxsubdiv`, try \
        increasing the keyword argument `maxsubdiv`."
    end

    return I, E
end

"""
    allocate_buffer(
        fct::FCT,
        domain::DOM;
        rule::EmbeddedCubature{D}=default_rule(domain),
        norm=LinearAlgebra.norm,
    ) where {D,FCT,DOM<:AbstractDomain{D}}

Allocate and return a heap buffer compatible with [`integrate`](@ref).

Passing this buffer through the `buffer` keyword can reduce memory allocations when
`integrate` is called repeatedly with compatible domain and value types.
"""
function allocate_buffer(
    fct::FCT,
    domain::DOM;
    rule::EmbeddedCubature{D}=default_rule(domain),
    norm=LinearAlgebra.norm,
) where {D,FCT,DOM<:AbstractDomain{D}}
    # Determine the type of elements returned by the embedded cubature.
    I, E = rule(fct, domain, norm)

    # Create a binary heap to store elements of the form (domain, I, E), where:
    # - `domain` is the current subdomain being processed.
    # - `I` is the integral value over the subdomain.
    # - `E` is the error estimate over the subdomain.
    # The heap is ordered by the maximum error (E) in descending order.
    buffer = BinaryHeap{Tuple{DOM,typeof(I),typeof(E)}}(
        Base.Order.By(last, Base.Order.Reverse)
    )

    return buffer
end

# paranoia about accumulated round-off, see re-sum from QuadGK.jl
"""
    resum(buffer; norm=LinearAlgebra.norm)

Re-sum integral and error contributions stored in `buffer`.

This is more expensive than a plain sum, but it uses the Kahan-Babuška-Neumaier [1]
summation algorithm to reduce floating-point round-off error.

[1] Klein, A. A Generalized Kahan-Babuška-Summation-Algorithm. Computing 76, 279-293 (2006).
https://doi.org/10.1007/s00607-005-0139-x
"""
function resum(
    buffer::BinaryHeap{Tuple{DOM,IT,ET}}; norm=LinearAlgebra.norm
) where {DOM,IT,ET}
    I = cᵢ = zero(IT)
    E = cₑ = zero(ET)
    for (_, I_dom, E_dom) in buffer.valtree
        tᵢ = I + I_dom
        if norm(I) ≥ norm(I_dom)
            cᵢ += (I - tᵢ) + I_dom
        else
            cᵢ += (I_dom - tᵢ) + I
        end
        I = tᵢ

        tₑ = E + E_dom
        if norm(E) ≥ norm(E_dom)
            cₑ += (E - tₑ) + E_dom
        else
            cₑ += (E_dom - tₑ) + E
        end
        E = tₑ
    end
    return I + cᵢ, E + cₑ
end
