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
    heap = allocate_buffer(f, s),
) where {N,T}
    return _integrate(f, s, quad, atol, rtol, maxsplit, norm, heap)
end

function _integrate(
    f,
    s::Simplex{N},
    quad::EmbeddedQuadrature{N},
    atol,
    rtol,
    maxsplit,
    norm,
    heap,
) where {N}
    nsplit = 0
    I, E = _integrate_with_error(f, s, quad)
    # a quick check to see if splitting is really needed
    if E < atol || E < rtol * norm(I) || nsplit >= maxsplit
        return I, E
    end
    # split is needed, so push the element to the heap and begin
    empty!(heap.valtree)
    push!(heap, (s, I, E))
    while E > atol && E > rtol * norm(I) && nsplit < maxsplit
        sc, Ic, Ec = pop!(heap)
        I -= Ic
        E -= Ec
        for child in subdivide(sc)
            # since the jacobian is constant, factor it out of the integration
            Inew, Enew = _integrate_with_error(f, child, quad)
            I += Inew
            E += Enew
            push!(heap, (child, Inew, Enew))
        end
        nsplit += 1
    end
    nsplit >= maxsplit && @warn "maximum number of steps reached"
    return I, E
end

function allocate_buffer(f, ::Simplex{N,T}) where {N,T}
    # try to infer the type of the element that will be returned by f
    S = Base.promote_op((x, w) -> f(x) * w, SVector{N,T}, T)
    isbitstype(S) || (@warn "non bitstype detected")
    # the heap of adaptive quadratures have elements of the form (s,I,E), where
    # I and E are the value and error estimate over the simplex s. The ordering
    # used is based the maximum error
    ord  = Base.Order.By(el -> -el[3])
    heap = BinaryHeap{Tuple{Simplex{N,T},S,T}}(ord)
    return heap
end
