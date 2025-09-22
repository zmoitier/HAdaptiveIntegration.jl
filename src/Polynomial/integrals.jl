function integral_monomials(::Type{<:Segment}, deg_tot_max::Int)
    return [[(k,) => 1//(k + 1)] for k in 0:deg_tot_max]
end

function integral_chebyshev(::Type{<:Segment}, deg_tot_max::Int)
    return [[(k,) => iseven(k) ? 1//(1 - k * k) : 0//1] for k in 0:deg_tot_max]
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
