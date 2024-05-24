function check_order(quad::Quadrature{N,T}, order::Int) where {N,T}
    exp_val = _exponent_to_value(N, order + 1)

    all_true = true
    for k in 0:order
        for α in _tot_deg(k, N)
            fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
            val = exp_val[α]
            if val == 0
                atol = √eps(T)
            else
                atol = 0
            end
            all_true = all_true & isapprox(quad(fct), T(val); atol = atol)
        end
    end

    all_true_next_order = true
    for α in _tot_deg(order + 1, N)
        fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
        val = exp_val[α]
        if val == 0
            atol = √eps(T)
        else
            atol = 0
        end
        all_true_next_order = all_true_next_order & isapprox(quad(fct), T(0); atol = atol)
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
    exponent_to_value = DefaultDict{Tuple{Int},Rational}(0 // 1)
    for k in 0:2:order
        exponent_to_value[(k,)] = 2 // (k + 1)
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
