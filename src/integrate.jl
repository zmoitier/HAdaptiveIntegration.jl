"""
    integrate(
        fct,
        domain::Domain{D,T};
        embedded_cubature::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain),
        subdiv_algo=default_subdivision(domain),
        buffer=nothing,
        norm=x -> LinearAlgebra.norm(x, Inf),
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),
        maxsubdiv=512 * 2^D,
    ) where {H,L,D,T}

Return `I` and `E` where `I` is the integral of the function `fct` over `domain` and `E` is
an error estimate.

## Arguments
- `fct`: a function that must take a `SVector{D,T}` to a return type `K`, with `K` must
   support the multiplication by a scalar of type `T` and the addition.
- `domain::Domain{D,T}`: the integration domain. Currently, we support [`segment`](@ref),
  [`triangle`](@ref), [`rectangle`](@ref), [`tetrahedron`](@ref), [`cuboid`](@ref), and
  d-dimensional [`simplex`](@ref).

## Optional arguments
- `ec::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain)`: the embedded cubature,
   each supported domain has a [`default_embedded_cubature`](@ref).
- `subdiv_algo=default_subdivision(domain)`: the subdivision algorithm, each domain has a
  [`default_subdivision`](@ref).
- `buffer=nothing`: heap use to do the adaptive algorithm, can be allocated using
   [`allocate_buffer`](@ref), which might result in performance gain if multiple call to
   integrate is perform.
- `norm=x -> LinearAlgebra.norm(x, Inf)`: norm used to estimate the error.
- `atol=zero(T)`: absolute tolerance.
- `rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T))`: relative tolerance.
- `maxsubdiv=512 * 2^D`: maximum number of subdivision.
"""
function integrate(
    fct,
    domain::Domain{D,T};
    embedded_cubature::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain),
    subdiv_algo=default_subdivision(domain),
    buffer=nothing,
    norm=x -> LinearAlgebra.norm(x, Inf),
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),
    maxsubdiv=512 * 2^D,
) where {H,L,D,T}
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
    nbsubdiv = 0
    I, E = ec(fct, domain, norm)

    # a quick check to see if splitting is really needed
    if (E < atol) || (E < rtol * norm(I)) || (nbsubdiv ≥ maxsubdiv)
        (nbsubdiv ≥ maxsubdiv) && @warn "maximum number of subdivide reached"
        return I, E
    end

    # split is needed, so initialize the heap
    heap = if isnothing(buffer)
        ord = Base.Order.By(el -> -el[3])
        BinaryHeap{Tuple{typeof(domain),typeof(I),typeof(E)}}(ord)
    else
        empty!(buffer.valtree)
        buffer
    end

    push!(heap, (domain, I, E))
    while (E > atol) && (E > rtol * norm(I)) && (nbsubdiv < maxsubdiv)
        sc, Ic, Ec = pop!(heap)
        I -= Ic
        E -= Ec
        for child in subdiv_algo(sc)
            Inew, Enew = ec(fct, child, norm)
            I += Inew
            E += Enew
            push!(heap, (child, Inew, Enew))
        end
        nbsubdiv += 1
    end

    (nbsubdiv ≥ maxsubdiv) && @warn "maximum number of subdivide reached"

    return I, E
end

"""
    allocate_buffer(
        fct, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)
    ) where {DOM<:Domain}

TBW
"""
function allocate_buffer(
    fct, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)
) where {DOM<:Domain}
    # Type of element that will be returned by the cubature. Pay the cost of single call to
    # figure this out.
    I, E = ec(fct, domain)
    # The heap of adaptive cubature have elements of the form (domain,I,E), where I and E
    # are the value and error estimate over `domain`. The ordering used is based on the
    # maximum error.
    heap = BinaryHeap{Tuple{DOM,typeof(I),typeof(E)}}(Base.Order.By(el -> -el[3]))
    return heap
end
