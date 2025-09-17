struct Errors{T}
    value::T
    absolute::T
    relative::T

    function Errors(value_num::T, value_ref::S) where {T,S}
        absolute = abs(value_num - value_ref)
        relative = absolute / abs(value_ref)
        return new{float(T)}(value_ref, absolute, relative)
    end
end

function Base.show(io::IO, e::Errors)
    @printf(io, "{v≈%.2e, ae≈%.1e, re≈%.1e}", e.value, e.absolute, e.relative)
end

function compute_error_monomials(
    tec::TabulatedEmbeddedCubature{DOM}, ::Type{T}; next_deg::Bool=false
) where {DOM,T}
    ec = embedded_cubature(tec, T)
    if next_deg
        return compute_error_monomials(ec, DOM, tec.order_high + 1, tec.order_low + 1)
    end
    return compute_error_monomials(ec, DOM, tec.order_high, tec.order_low)
end

function compute_error_monomials(
    ec::EmbeddedCubature{D,T}, ::Type{DOM}, order_high::Int, order_low::Int
) where {D,T,DOM}
    exponent2values = integral_monomials(DOM, order_high)

    α2e_h = Dict{Int,Dict{NTuple{D,Int},Errors{T}}}()
    α2e_l = Dict{Int,Dict{NTuple{D,Int},Errors{T}}}()

    L = length(ec.weights_low)
    for k in 0:order_low
        tmp_hi = Dict{NTuple{D,Int},Errors{T}}()
        tmp_lo = Dict{NTuple{D,Int},Errors{T}}()
        for (e, v) in exponent2values[k + 1]
            vh, vl = zero(T), zero(T)
            for (i, node) in enumerate(ec.nodes)
                y = prod(node .^ e)
                vh += ec.weights_high[i] * y
                if i ≤ L
                    vl += ec.weights_low[i] * y
                end
            end
            tmp_hi[e] = Errors(vh, v)
            tmp_lo[e] = Errors(vl, v)
        end
        α2e_h[k] = tmp_hi
        α2e_l[k] = tmp_lo
    end
    for k in (order_low + 1):order_high
        tmp_hi = Dict{NTuple{D,Int},Errors{T}}()
        for (e, v) in exponent2values[k + 1]
            vh = zero(T)
            for (i, node) in enumerate(ec.nodes)
                vh += ec.weights_high[i] * prod(node .^ e)
            end
            tmp_hi[e] = Errors(vh, v)
        end
        α2e_h[k] = tmp_hi
    end

    return α2e_h, α2e_l
end

function integral_monomials(::Type{<:Segment}, deg_tot_max::Int)
    return [[(i,) => 1//(i + 1)] for i in 0:deg_tot_max]
end

function integral_monomials(::Type{<:Simplex{D}}, deg_tot_max::Int) where {D}
    @assert (D > 0) && (deg_tot_max ≥ 0) "must have `D > 0` and `k_max ≥ 0`."

    exponent2values = [[(i,) => 1//prod((i + 1):(i + D))] for i in 0:deg_tot_max]

    for d in 2:D
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, α2v) in zip(countfrom(0), exponent2values)
            for (α, v) in α2v
                push!(new[k + 1], (0, α...) => v)
                t = 1
                for n in 1:(deg_tot_max - k)
                    t *= n//(n + k + D)
                    push!(new[k + 1 + n], (n, α...) => t * v)
                end
            end
        end
        exponent2values = new
    end

    return exponent2values
end

function integral_monomials(::Type{<:Orthotope{D}}, deg_tot_max::Int) where {D}
    @assert (D > 0) && (deg_tot_max ≥ 0) "must have `dim > 0` and `k_max ≥ 0`."

    exponent2values = [[(i,) => 1//(i + 1)] for i in 0:deg_tot_max]

    for d in 2:D
        new = [Vector{Pair{NTuple{d,Int},Rational{Int}}}() for _ in 0:deg_tot_max]
        for (k, α2v) in zip(countfrom(0), exponent2values)
            for (α, v) in α2v
                for n in 0:(deg_tot_max - k)
                    push!(new[k + 1 + n], (n, α...) => v//(n + 1))
                end
            end
        end
        exponent2values = new
    end

    return exponent2values
end
