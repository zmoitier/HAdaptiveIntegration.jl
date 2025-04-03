"""
    validate_orders(
        ec::EmbeddedCubature{D,T},
        ::Type{DOM},
        order_high::Int,
        order_low::Int;
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {D,T,DOM<:AbstractDomain}

Return `true` if the embedded cubature on the **reference** domain of `DOM` integrate
exactly (within tolerance) the monomials up to degree `order_high` for the high order
cubature and `order_low` for the low order cubature. Else return `false`.
"""
function validate_orders(
    ec::EmbeddedCubature{D,T},
    ::Type{DOM},
    order_high::Int,
    order_low::Int;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T,DOM<:AbstractDomain}
    if DOM <: Orthotope
        val_ref = integral_monomial_orthotope(D, order_high)
    elseif DOM <: Simplex
        val_ref = integral_monomial_simplex(D, order_high)
    else
        @error "Unsupported domain type: $DOM)"
        return false
    end

    return _validate_orders(ec, val_ref, order_high, order_low; atol=atol, rtol=rtol)
end

function _validate_orders(
    ec::EmbeddedCubature{D,T},
    val_ref::Vector{Vector{Pair{NTuple{D,Int},Rational{Int}}}},
    order_high::Int,
    order_low::Int;
    atol::T,
    rtol::T,
) where {D,T}
    for k in 0:order_high
        for (idx, val_ref) in val_ref[k + 1]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, idx))
            Iₕ, Iₗ = ec(fct)

            vᵣ = T(val_ref)
            if k ≤ order_low
                if !isapprox(Iₗ, vᵣ; atol=atol, rtol=rtol)
                    @error "fail to integrate the low order cubature within tolerance at degree = $idx. Computed value: $Iₗ, Reference value: $vᵣ."
                    return false
                end
            end

            if !isapprox(Iₕ, vᵣ; atol=atol, rtol=rtol)
                @error "fail to integrate the high order cubature within tolerance at degree = $idx. Computed value: $Iₗ, Reference value: $vᵣ."
                return false
            end
        end
    end

    return true
end

function (ec::EmbeddedCubature{D,T})(fct) where {D,T}
    H, L = length(ec.weights_high), length(ec.weights_low)

    v = fct(ec.nodes[1])
    Iₗₒ = ec.weights_low[1] * v
    Iₕᵢ = ec.weights_high[1] * v
    for i in 2:L
        v = fct(ec.nodes[i])
        Iₗₒ += ec.weights_low[i] * v
        Iₕᵢ += ec.weights_high[i] * v
    end
    for i in (L + 1):H
        Iₕᵢ += ec.weights_high[i] * fct(ec.nodes[i])
    end
    return Iₕᵢ, Iₗₒ
end

"""
    integral_monomial_orthotope(dim::Int, k_max::Int)

Return the values of monomial's integral over the **reference** orthotope. It return a
`Vector{ Vector{ Pair{ Ntuple{dim,Int}, Rational{Int} } } }`. The outer vector is index by
the `total degree + 1`, for the total degree form 0 to `tot_deg_max`. The inner vector
contain `Pair{ Ntuple{dim,Int}, Rational{Int} }` where the `Ntuple{D,Int}` is the
multi-index of the monomial and `Rational{Int}` is the value of the integral.
"""
function integral_monomial_orthotope(dim::Int, k_max::Int)
    @assert (dim > 0) && (k_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    mlt_idx_val = [[(i,) => 1//(i + 1)] for i in 0:k_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:k_max]
        for (k, mi_val) in zip(Iterators.countfrom(0), mlt_idx_val)
            for (mi, val) in mi_val
                for n in 0:(k_max - k)
                    push!(new[k + 1 + n], (n, mi...) => val//(n + 1))
                end
            end
        end
        mlt_idx_val = new
    end

    return mlt_idx_val
end

"""
    integral_monomial_simplex(dim::Int, k_max::Int)

Return the values of monomial's integral over the **reference** simplex. It return a
`Vector{ Vector{ Pair{ Ntuple{dim,Int}, Rational{Int} } } }`. The outer vector is index by
the `total degree + 1`, for the total degree form 0 to `tot_deg_max`. The inner vector
contain `Pair{ Ntuple{dim,Int}, Rational{Int} }` where the `Ntuple{D,Int}` is the
multi-index of the monomial and `Rational{Int}` is the value of the integral.
"""
function integral_monomial_simplex(dim::Int, k_max::Int)
    @assert (dim > 0) && (k_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    mlt_idx_val = [[(i,) => 1//prod((i + 1):(i + dim))] for i in 0:k_max]

    for d in 2:dim
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:k_max]
        for (k, mi_val) in zip(Iterators.countfrom(0), mlt_idx_val)
            for (mi, val) in mi_val
                push!(new[k + 1], (0, mi...) => val)
                v = 1
                for n in 1:(k_max - k)
                    v *= n//(n + k + dim)
                    push!(new[k + 1 + n], (n, mi...) => val * v)
                end
            end
        end
        mlt_idx_val = new
    end

    return mlt_idx_val
end
