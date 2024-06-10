"""
    check_order(quad::Quadrature{N,T}, order::Int)::Bool where {N,T}

Check if the quadrature rule `quad` has order `order` by integrating monomials.
"""
function check_order(quad::Quadrature{N,T}, order::Int)::Bool where {N,T}
    rtol = 10 * eps(T)
    exp_val = _exponent_to_value(order + 1, N)

    all_true = true
    for k in 0:order
        for α in _multi_index(k, N)
            fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
            is_integrated = isapprox(quad(fct), T(exp_val[α]); rtol = rtol)
            all_true = all_true & is_integrated

            @debug α is_integrated abs(quad(fct) / T(exp_val[α]) - T(1))
        end
    end

    all_true_next_order = true
    for α in _multi_index(order + 1, N)
        fct = x -> prod(xᵢ^αᵢ for (xᵢ, αᵢ) in zip(x, α))
        is_integrated = isapprox(quad(fct), T(exp_val[α]); rtol = rtol)
        all_true_next_order = all_true_next_order & is_integrated

        @debug α is_integrated abs(quad(fct) / T(exp_val[α]) - T(1))
    end

    return all_true & !all_true_next_order
end

"""
    _exponent_to_value(tot_deg_max::Int, N::Int)::Dict{NTuple{N,Int64},Rational{Int64}}

Generate a dictionary mapping multi-indices to the value of the integral of the
monomial on the reference simplex.
"""
function _exponent_to_value(tot_deg_max::Int, N::Int)::Dict{NTuple{N,Int64},Rational{Int64}}
    exponent_to_value =
        Dict{NTuple{N,Int64},Rational{Int64}}(Tuple(0 for _ in 1:N) => 1 // factorial(N))

    for k in 1:tot_deg_max
        for exponent in _multi_index(k, N)
            decr_exp, i = _decr_multi_idx(exponent)
            exponent_to_value[exponent] =
                exponent_to_value[decr_exp] * Rational(exponent[i], sum(exponent) + N)
        end
    end

    return exponent_to_value
end

"""
    _decr_multi_idx(exponent::NTuple{N,Int64})::Tuple{NTuple{N,Int64},Int} where {N}

Decrement the first positive entry of the multi-index `exponent` and return the new
multi-index and the position of the first positive entry.
"""
function _decr_multi_idx(exponent::NTuple{N,Int64})::Tuple{NTuple{N,Int64},Int} where {N}
    i = findfirst(x -> x > 0, exponent)
    @assert i !== nothing "No positive entry in the exponent"

    decr_exp = Tuple(j ≠ i ? e : e - 1 for (j, e) in enumerate(exponent))
    return decr_exp, i
end

"""
    _multi_index(tot_deg::Int, N::Int)::Vector{NTuple{N,Int}}

Generate all multi-indices of length `N` with sum of entries equal to `tot_deg`.
"""
function _multi_index(tot_deg::Int, N::Int)::Vector{NTuple{N,Int}}
    if N ≤ 0
        return []
    end

    if N == 1
        return [(tot_deg,)]
    end

    return [(d, t...) for d in tot_deg:-1:0 for t in _multi_index(tot_deg - d, N - 1)]
end
