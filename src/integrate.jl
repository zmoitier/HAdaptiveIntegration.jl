"""
    integrate(f, s::Simplex[, quad; atol = 0, rtol = 0, maxsplit, norm])
"""
function integrate(
    f,
    s::Simplex{N,T},
    quad::EmbeddedQuadrature{N,T} = default_quadrature(s);
    atol = zero(T),
    rtol = atol == zero(T) ? sqrt(eps(T)) : zero(T),
    maxsplit = 1000,
    norm = LinearAlgebra.norm,
    heap = nothing,
) where {N,T}
    return _integrate(f, s, quad, atol, rtol, maxsplit, norm, heap)
end

function _integrate(
    f,
    s::Simplex{N,T},
    quad::EmbeddedQuadrature{N,T},
    atol,
    rtol,
    maxsplit,
    norm,
    buf,
) where {N,T}
    nsplit = 0
    I, E = quad(f, s)
    # a quick check to see if splitting is really needed
    if E < atol || E < rtol * norm(I) || nsplit >= maxsplit
        return I, E
    end
    # split is needed, so prepare heap if needed, push the element to the heap
    # and begin
    heap = if isnothing(buf)
        ord = Base.Order.By(el -> -el[3])
        BinaryHeap{Tuple{Simplex{N,T},typeof(I),typeof(E)}}(ord)
    else
        empty!(buf.valtree)
        buf
    end
    push!(heap, (s, I, E))
    while E > atol && E > rtol * norm(I) && nsplit < maxsplit
        sc, Ic, Ec = pop!(heap)
        I -= Ic
        E -= Ec
        for child in subdivide(sc)
            # since the jacobian is constant, factor it out of the integration
            Inew, Enew = quad(f, child)
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
    f,
    s::Simplex{N,T},
    quad::EmbeddedQuadrature{N,T} = default_quadrature(s),
) where {N,T}
    # type of element that will be returned by by quad. Pay the cost of single
    # call to figure this out
    I, E = _integrate_with_error(f, s, quad)
    # the heap of adaptive quadratures have elements of the form (s,I,E), where
    # I and E are the value and error estimate over the simplex s. The ordering
    # used is based the maximum error
    ord  = Base.Order.By(el -> -el[3])
    heap = BinaryHeap{Tuple{Simplex{N,T},typeof(I),typeof(E)}}(ord)
    return heap
end
