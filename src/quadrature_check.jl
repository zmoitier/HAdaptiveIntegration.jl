function check_order(quad::Quadrature{N,T}, order::Int) where {N,T}
    rtol = 10 * eps(T)
    exp_val = _exponent_to_value(N, order + 1)

    all_true = true
    for k in 0:order
        for α in _tot_deg(k, N)
            fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
            is_integrated = isapprox(quad(fct), T(exp_val[α]); rtol = rtol)
            all_true = all_true & is_integrated

            @debug α is_integrated abs(quad(fct) / T(exp_val[α]) - T(1))
        end
    end

    all_true_next_order = true
    for α in _tot_deg(order + 1, N)
        fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
        is_integrated = isapprox(quad(fct), T(exp_val[α]); rtol = rtol)
        all_true_next_order = all_true_next_order & is_integrated

        @debug α is_integrated abs(quad(fct) / T(exp_val[α]) - T(1))
    end

    return all_true & !all_true_next_order
end

function _exponent_to_value(dim::Int, order::Int)
    if dim == 1
        return _seg_ref_val(order)
    else
        error("Unknown reference values in dimension $dim")
    end
end

function _seg_ref_val(order::Int)
    exponent_to_value = Dict{Tuple{Int},Rational{Int64}}()
    for k in 0:order
        exponent_to_value[(k,)] = 1 // (k + 1)
    end
    return exponent_to_value
end

function _tot_deg(tot_deg::Int, space_dim::Int)
    if space_dim ≤ 0
        return []
    end

    if space_dim == 1
        return [(tot_deg,)]
    end

    return [(d, t...) for d in tot_deg:-1:0 for t in _tot_deg(tot_deg - d, space_dim - 1)]
end
