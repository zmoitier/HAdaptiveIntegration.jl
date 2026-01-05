function integral_monomials_rule(
    tec::TabulatedEmbeddedCubature{DOM}, ::Type{T}; next_deg::Bool=false
) where {DOM,T}
    ec = embedded_cubature(tec, T)
    if next_deg
        return integral_monomials_rule(ec, tec.order_high + 1, tec.order_low + 1)
    end
    return integral_monomials_rule(ec, tec.order_high, tec.order_low)
end

function integral_monomials_rule(
    ec::EmbeddedCubature{D,T}, order_high::Int, order_low::Int
) where {D,T}
    exponents = multi_index(D, order_high)

    α2v_h = Vector{Vector{Pair{NTuple{D,Int},T}}}()
    α2v_l = Vector{Vector{Pair{NTuple{D,Int},T}}}()

    L = length(ec.weights_low)
    for k in 0:order_low
        tmp_hi = Vector{Pair{NTuple{D,Int},T}}()
        tmp_lo = Vector{Pair{NTuple{D,Int},T}}()
        for e in exponents[k + 1]
            vh, vl = zero(T), zero(T)
            for (i, node) in enumerate(ec.nodes)
                y = prod(node .^ e)
                vh += ec.weights_high[i] * y
                if i ≤ L
                    vl += ec.weights_low[i] * y
                end
            end
            push!(tmp_hi, e => vh)
            push!(tmp_lo, e => vl)
        end
        push!(α2v_h, tmp_hi)
        push!(α2v_l, tmp_lo)
    end
    for k in (order_low + 1):order_high
        tmp_hi = Vector{Pair{NTuple{D,Int},T}}()
        for e in exponents[k + 1]
            vh = zero(T)
            for (i, node) in enumerate(ec.nodes)
                vh += ec.weights_high[i] * prod(node .^ e)
            end
            push!(tmp_hi, e => vh)
        end
        push!(α2v_h, tmp_hi)
    end

    return α2v_h, α2v_l
end

function multi_index(dim::Int, deg_tot_max::Int)
    @assert (dim > 0) && (deg_tot_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    multi_idx = [[(i,)] for i in 0:deg_tot_max]
    for d in 2:dim
        tmp = [Vector{NTuple{d,Int}}() for _ in 0:deg_tot_max]
        for (k, αs) in zip(countfrom(0), multi_idx)
            for α in αs
                for n in 0:(deg_tot_max - k)
                    push!(tmp[k + 1 + n], (n, α...))
                end
            end
        end
        multi_idx = tmp
    end

    return multi_idx
end
