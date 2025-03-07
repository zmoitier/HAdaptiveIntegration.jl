"""
    integrate(
        fct,
        domain::Domain{D,T},
        ec::EmbeddedCubature{H,L,D,T}=default_embedded_cubature(domain);
        subdiv_algo=default_subdivision(domain),
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : sqrt(eps(T)),
        maxsplit=D * 1024,
        norm=LinearAlgebra.norm,
        buffer=nothing,
    ) where {H,L,D,T<:Real}

TBW
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
    maxsplit=D * 1024,
) where {H,L,D,T<:Real}
    return _integrate(
        fct, domain, embedded_cubature, subdiv_algo, buffer, norm, atol, rtol, maxsplit
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
    maxsplit,
) where {FCT,DOM}
    nsplit = 0
    I, E = ec(fct, domain, norm)

    # a quick check to see if splitting is really needed
    if (E < atol) || (E < rtol * norm(I)) || (nsplit >= maxsplit)
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
    while (E > atol) && (E > rtol * norm(I)) && (nsplit < maxsplit)
        sc, Ic, Ec = pop!(heap)
        I -= Ic
        E -= Ec
        for child in subdiv_algo(sc)
            Inew, Enew = ec(fct, child, norm)
            I += Inew
            E += Enew
            push!(heap, (child, Inew, Enew))
        end
        nsplit += 1
    end

    (nsplit >= maxsplit) && @warn "maximum number of subdivide reached"

    return I, E
end

"""
    allocate_buffer(
        fct::Function, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)
    ) where {DOM<:Domain}

TBW
"""
function allocate_buffer(
    fct::Function, domain::DOM, ec::EmbeddedCubature=default_embedded_cubature(domain)
) where {DOM<:Domain}
    # type of element that will be returned by quad. Pay the cost of single
    # call to figure this out
    I, E = ec(fct, domain)
    # the heap of adaptive quadratures have elements of the form (s,I,E), where
    # I and E are the value and error estimate over the simplex s. The ordering
    # used is based on the maximum error
    heap = BinaryHeap{Tuple{DOM,typeof(I),typeof(E)}}(Base.Order.By(el -> -el[3]))
    return heap
end
