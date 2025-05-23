"""
    integrate(
        fct,
        domain::AbstractDomain{D,T};
        embedded_cubature::EmbeddedCubature{D,T}=default_embedded_cubature(domain),
        subdiv_algo=default_subdivision(domain),
        buffer=nothing,
        norm=x -> LinearAlgebra.norm(x, Inf),
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),
        maxsubdiv=8192 * 2^D,
    ) where {D,T}

Return `I` and `E` where `I` is the integral of the function `fct` over `domain` and `E` is
an error estimate.

## Arguments
- `fct`: a function that must take a `SVector{D,T}` to a return type `K`, with `K` must
   support the multiplication by a scalar of type `T` and the addition.
- `domain::AbstractDomain{D,T}`: the integration domain. Currently, we support
  [`Triangle`](@ref), [`Rectangle`](@ref), [`Tetrahedron`](@ref), [`Cuboid`](@ref),
  `d`-dimensional [`Simplex`](@ref), and `d`-dimensional [`Orthotope`](@ref).

## Optional arguments
- `embedded_cubature::EmbeddedCubature{D,T}=default_embedded_cubature(domain)`: the embedded cubature,
   each supported domain has a [`default_embedded_cubature`](@ref).
- `subdiv_algo=default_subdivision(domain)`: the subdivision algorithm, each domain has a
  [`default_subdivision`](@ref).
- `buffer=nothing`: heap use to do the adaptive algorithm, can be allocated using
   [`allocate_buffer`](@ref), which might result in performance gain if multiple call to
   integrate is perform.
- `norm=x -> LinearAlgebra.norm(x, Inf)`: norm used to estimate the error.
- `atol=zero(T)`: absolute tolerance.
- `rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T))`: relative tolerance.
- `maxsubdiv=8192 * 2^D`: maximum number of subdivision.
"""
function integrate(
    fct,
    domain::AbstractDomain{D,T};
    embedded_cubature::EmbeddedCubature{D,T}=default_embedded_cubature(domain),
    subdiv_algo=default_subdivision(domain),
    buffer=nothing,
    norm=x -> norm(x, Inf),
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),
    maxsubdiv=8192 * 2^D,
) where {D,T}
    return _integrate(
        fct, domain, embedded_cubature, subdiv_algo, buffer, norm, atol, rtol, maxsubdiv
    )
end

function _integrate(
    fct::FCT,
    domain::DOM,
    ec::EmbeddedCubature,
    subdiv_algo,
    buffer,
    norm,
    atol,
    rtol,
    maxsubdiv,
) where {FCT,DOM}
    nb_subdiv = 0
    I, E = ec(fct, domain, norm)

    # a quick check to see if splitting is really needed
    if (E < atol) || (E < rtol * norm(I)) || (nb_subdiv ≥ maxsubdiv)
        return I, E
    end

    # split is needed, so initialize the heap
    heap = if isnothing(buffer)
        BinaryHeap{Tuple{typeof(domain),typeof(I),typeof(E)}}(
            Base.Order.By(last, Base.Order.Reverse)
        )
    else
        empty!(buffer.valtree)
        buffer
    end

    push!(heap, (domain, I, E))
    while (E > atol) && (E > rtol * norm(I)) && (nb_subdiv < maxsubdiv)
        domain, I_dom, E_dom = pop!(heap)
        I -= I_dom
        E -= E_dom
        for child in subdiv_algo(domain)
            I_child, E_child = ec(fct, child, norm)
            I += I_child
            E += E_child
            push!(heap, (child, I_child, E_child))
        end
        nb_subdiv += 1
    end

    if nb_subdiv ≥ maxsubdiv
        @warn "maximum number of subdivide reached, try increasing the keyword argument `maxsubdiv=$maxsubdiv`."
    end

    return I, E
end

"""
    allocate_buffer(fct, domain, ec=default_embedded_cubature(domain))

Allocate and return a buffer that can be passed to the [`integrate`](@ref) function
to improve performance by reducing memory allocations when `integrate` is called
multiple times.
"""
function allocate_buffer(
    fct, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)
) where {DOM<:AbstractDomain}
    # Determine the type of elements returned by the embedded cubature.
    I, E = ec(fct, domain)

    # Create a binary heap to store elements of the form (domain, I, E), where:
    # - `domain` is the current subdomain being processed.
    # - `I` is the integral value over the subdomain.
    # - `E` is the error estimate over the subdomain.
    # The heap is ordered by the maximum error (E) in descending order.
    heap = BinaryHeap{Tuple{typeof(domain),typeof(I),typeof(E)}}(
        Base.Order.By(last, Base.Order.Reverse)
    )

    return heap
end
