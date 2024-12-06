"""
    integrate(f, s::Simplex[, quad, subdiv_algo; atol = 0, rtol = 0, maxsplit,
    norm])
"""
function integrate(
    fct,
    domain::Domain{N,T},
    quad::EmbeddedQuadrature{N,T}=default_quadrature(domain),
    subdiv_algo=default_subdivision(domain);
    atol=zero(T),
    rtol=atol == zero(T) ? sqrt(eps(T)) : zero(T),
    maxsplit=typemax(Int),
    norm=LinearAlgebra.norm,
    heap=nothing,
) where {N,T}
    return _integrate(fct, domain, quad, subdiv_algo, atol, rtol, maxsplit, norm, heap)
end

function _integrate(
    fct::F,
    domain,
    quad::EmbeddedQuadrature,
    subdiv_algo,
    atol,
    rtol,
    maxsplit,
    norm,
    buffer,
) where {F}
    nsplit = 0
    I, E = quad(fct, domain)

    # a quick check to see if splitting is really needed
    if E < atol || E < rtol * norm(I) || nsplit >= maxsplit
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
    while E > atol && E > rtol * norm(I) && nsplit < maxsplit
        sc, Ic, Ec = pop!(heap)
        I -= Ic
        E -= Ec
        for child in subdiv_algo(sc)
            Inew, Enew = quad(fct, child)
            I += Inew
            E += Enew
            push!(heap, (child, Inew, Enew))
        end
        nsplit += 1
    end

    nsplit >= maxsplit && @warn "maximum number of steps reached"

    return I, E
end

function allocate_buffer(
    fct, domain::D, quad::EmbeddedQuadrature=default_quadrature(domain)
) where {D}
    # type of element that will be returned by quad. Pay the cost of single
    # call to figure this out
    I, E = quad(fct, domain)
    # the heap of adaptive quadratures have elements of the form (s,I,E), where
    # I and E are the value and error estimate over the simplex s. The ordering
    # used is based on the maximum error
    ord = Base.Order.By(el -> -el[3])
    heap = BinaryHeap{Tuple{D,typeof(I),typeof(E)}}(ord)
    return heap
end
