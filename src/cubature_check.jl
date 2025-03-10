"""
    check_order(
        tec::TabulatedEmbeddedCubature,
        domain::Domain{D,T};
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {D,T}

Return 0 if the tabulated embedded cubature on the **reference** domain of `domain`
integrate exactly (within tolerance) the monomials up to degree `order_high` for the high
oder cubature and `order_low` for the low order cubature. Else return 1.
"""
function check_order(
    tec::TabulatedEmbeddedCubature,
    domain::Domain{D,T};
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {D,T}
    return check_order(
        embedded_cubature(tec, T),
        domain,
        tec.order_high,
        tec.order_low;
        atol=atol,
        rtol=rtol,
    )
end

"""
    check_order(
        ec::EmbeddedCubature{H,L,D,T},
        domain::Domain{D,T},
        order_high::Int,
        order_low::Int;
        atol=zero(T),
        rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
    ) where {H,L,D,T}

Return 0 if the embedded cubature on the **reference** domain of `domain` integrate exactly
(within tolerance) the monomials up to degree `order_high` for the high oder cubature and
`order_low` for the low order cubature. Else return 1.
"""
function check_order(
    ec::EmbeddedCubature{H,L,D,T},
    domain::Domain{D,T},
    order_high::Int,
    order_low::Int;
    atol=zero(T),
    rtol=(atol > zero(T)) ? zero(T) : 10 * eps(T),
) where {H,L,D,T}
    val_ref = integral_monomial(domain, order_high)

    val_lo::Vector{Vector{T}} = []
    val_hi::Vector{Vector{T}} = []
    for i in 1:(order_low + 1)
        tmp_lo = Vector{T}()
        tmp_hi = Vector{T}()
        for (idx, _) in val_ref[i]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, idx))
            Ih, Il = ec(fct)
            push!(tmp_lo, Il)
            push!(tmp_hi, Ih)
        end
        push!(val_lo, tmp_lo)
        push!(val_hi, tmp_hi)
    end
    for i in (order_low + 2):(order_high + 1)
        tmp_hi = Vector{T}()
        for (idx, _) in val_ref[i]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, idx))
            Ih, _ = ec(fct)
            push!(tmp_hi, Ih)
        end
        push!(val_hi, tmp_hi)
    end

    if _isapprox(val_ref, val_hi, atol, rtol) == 1
        @error "integration fail for high order cubature"
        return 1
    end
    if _isapprox(val_ref, val_lo, atol, rtol) == 1
        @error "integration fail for low order cubature"
        return 1
    end
    return 0
end

function (ec::EmbeddedCubature{H,L,D,T})(fct) where {H,L,D,T}
    v = fct(ec.nodes[1])
    I_low = v * ec.weights_low[1]
    I_high = v * ec.weights_high[1]
    for i in 2:L
        v = fct(ec.nodes[i])
        I_low += v * ec.weights_low[i]
        I_high += v * ec.weights_high[i]
    end
    for i in (L + 1):H
        I_high += fct(ec.nodes[i]) * ec.weights_high[i]
    end
    return I_high, I_low
end

function _isapprox(
    val_ref::Vector{Vector{Pair{NTuple{D,Int},Rational{Int}}}},
    val_num::Vector{Vector{T}},
    atol::T,
    rtol::T,
) where {D,T<:Real}
    for (i, (idx_ref, num)) in enumerate(zip(val_ref, val_num))
        for ((_, vr), vn) in zip(idx_ref, num)
            if !isapprox(vn, vr; atol=atol, rtol=rtol)
                @error "fail to integrate within tolerance at total degree = $(i-1)"
                return 1
            end
        end
    end
    @info "integrate within tolerance up to total degree = $(length(val_num) - 1)"
    return 0
end

"""
    integral_monomial(domain::Domain{D,T}, tot_deg_max::Int) where {D,T}

Return the values of monomial's integral over the **reference** domain. It return a
`Vector{ Vector{ Pair{ Ntuple{dim,Int}, Rational{Int} } } }`. The outer vector is index by
the `total degree + 1`, for the total degree form 0 to `tot_deg_max`. The inner vector
contain `Pair{ Ntuple{dim,Int}, Rational{Int} }` where the `Ntuple{D,Int}` is the
multi-index of the monomial and `Rational{Int}` is the value of the integral.
"""
function integral_monomial(::Orthotope{D,T}, total_degree_max::Int) where {D,T}
    return integral_monomial_orthotope(D, total_degree_max)
end

function integral_monomial(::Simplex{D,T,N}, total_degree_max::Int) where {D,T,N}
    return integral_monomial_simplex(D, total_degree_max)
end

function integral_monomial_orthotope(d::Int, tdm::Int)
    if d ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(i,) => 1//(i + 1)] for i in 0:tdm]

    if d == 1
        return indexes
    end

    for n in 2:d
        tmp::Vector{Vector{Pair{NTuple{n,Int},Rational{Int}}}} = [[] for _ in 0:tdm]
        for (i, vec) in enumerate(indexes)
            for (idx, val) in vec
                for k in 0:(tdm - i + 1)
                    push!(tmp[i + k], (k, idx...) => val//(k + 1))
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

function integral_monomial_simplex(d::Int, tdm::Int)
    if d ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(i,) => 1//prod((i + 1):(i + d))] for i in 0:tdm]

    if d == 1
        return indexes
    end

    for n in 2:d
        tmp::Vector{Vector{Pair{NTuple{n,Int},Rational{Int}}}} = [[] for _ in 0:tdm]
        for (i, vec) in enumerate(indexes)
            for (idx, val) in vec
                push!(tmp[i], (0, idx...) => val)
                v = 1
                for k in 1:(tdm - i + 1)
                    v *= k//(k + i - 1 + d)
                    push!(tmp[i + k], (k, idx...) => val * v)
                end
            end
        end
        indexes = tmp
    end

    return indexes
end
