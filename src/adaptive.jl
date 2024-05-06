Base.@kwdef struct AdaptiveQuadrature{N,T}
    qrule::EmbeddedQuadrature{N,T}
    atol::T = 0
    rtol::T = atol > 0 ? 0 : sqrt(eps(T))
    maxsplit::Int = typemax(Int)
end

function integrate_with_error(
    f,
    t::Simplex{N},
    q::AdaptiveQuadrature{N},
    heap = allocate_buffer(f, q),
) where {N}
    μ = measure(t)
    g = (s) -> f(t(s))
    I, E = integrate_with_error(f, q.qrule)
    nsplit = 0
    # a quick check to see if splitting is really needed
    if E < q.atol || E < q.rtol * norm(I) || nsplit >= q.maxsplit
        return I, E
    end
    # split is needed, so push the element to the heap and begin
    empty!(heap.valtree)
    T = eltype(heap).parameters[1]
    t = T(vertices(domain(q)))
    push!(heap, t => (I, E))
    I, E = _integrate_with_error!(f, heap, I, E, nsplit, q)
    return I, E
end

function _integrate_with_error!(f::F, heap, I, E, nsplit, q) where {F}
    while E > q.atol && E > q.rtol * norm(I) && nsplit < q.maxsplit
        t, (Ic, Ec) = pop!(heap)
        I -= Ic
        E -= Ec
        for child in split(t, q)
            # since the jacobian is constant, factor it out of the integration
            g = (s) -> f(child(s))
            μ = measure(child)
            Inew, Enew = μ .* integrate_with_error(g, q.qrule)
            I += Inew
            E += Enew
            push!(heap, child => (Inew, Enew))
        end
        nsplit += 1
    end
    nsplit >= q.maxsplit && @warn "maximum number of steps reached"
    return I, E
end

function allocate_buffer(f, q::AdaptiveQuadrature{N,T}) where {N,T}
    # try to infer the type of the element that will be returned by f
    S = return_type(f, eltype(q.nodes))
    TS = Base.promote_op(*, S, T)
    isbitstype(TS) || (@warn "non bitstype detected")
    D = Simplex{N,T,N+1}
    # the heap of adaptive quadratures have elements of the form s => (I,E),
    # where I and E are the value and error estimate over the simplex s. The
    # ordering used is based on the negative of the error
    ord = Base.Order.By(el -> -el[2][2])
    heap = BinaryHeap{Pair{D,Tuple{TS,T}}}(ord)
    return heap
end
