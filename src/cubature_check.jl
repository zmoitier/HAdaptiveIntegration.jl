function check_order(
    ecr::EmbeddedCubatureRaw,
    domain_type::DataType;
    atol::Union{Real,Nothing}=nothing,
    rtol::Union{Real,Nothing}=nothing,
)
    D = domain_type.parameters[1]
    T = domain_type.parameters[2]

    ol, oh = ecr.order_low, ecr.order_high
    if domain_type <: Simplex
        val_ref = simplex_int_val(D, oh)
    elseif domain_type <: Orthotope
        val_ref = orthotope_int_val(D, oh)
    else
        @error "unknown monomial's integrals on the reference domain of $domain_type."
    end

    ec = embedded_cubature_from_raw(ecr, T)

    val_lo::Vector{Vector{Pair{NTuple{D,Int},T}}} = []
    val_hi::Vector{Vector{Pair{NTuple{D,Int},T}}} = []
    for i in 1:(ol + 1)
        tmp_lo::Vector{Pair{NTuple{D,Int},T}} = []
        tmp_hi::Vector{Pair{NTuple{D,Int},T}} = []
        for (idx, _) in val_ref[i]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, idx))
            Ih, Il = ec(fct)
            push!(tmp_lo, idx => Il)
            push!(tmp_hi, idx => Ih)
        end
        push!(val_lo, tmp_lo)
        push!(val_hi, tmp_hi)
    end
    for i in (ol + 2):(oh + 1)
        tmp_hi::Vector{Pair{NTuple{D,Int},T}} = []
        for (idx, _) in val_ref[i]
            fct = x -> prod(xᵢ^eᵢ for (xᵢ, eᵢ) in zip(x, idx))
            Ih, _ = ec(fct)
            push!(tmp_hi, idx => Ih)
        end
        push!(val_hi, tmp_hi)
    end

    if isnothing(atol)
        atol = zero(T)
    end
    if isnothing(rtol)
        rtol = (atol > zero(T)) ? zero(T) : 10 * eps(T)
    end

    return nothing
end

function orthotope_int_val(dim::Int, tot_deg_max::Int)
    if dim ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(i,) => 1//(i + 1)] for i in 0:tot_deg_max]

    if dim == 1
        return indexes
    end

    for n in 2:dim
        tmp::Vector{Vector{Pair{NTuple{n,Int},Rational{Int}}}} = [[] for _ in 0:tot_deg_max]
        for (i, vec) in enumerate(indexes)
            for (idx, val) in vec
                for k in 0:(tot_deg_max - i + 1)
                    push!(tmp[i + k], (k, idx...) => val//(k + 1))
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

function simplex_int_val(dim::Int, tot_deg_max::Int)
    if dim ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(i,) => 1//prod((i + 1):(i + dim))] for i in 0:tot_deg_max]

    if dim == 1
        return indexes
    end

    for n in 2:dim
        tmp::Vector{Vector{Pair{NTuple{n,Int},Rational{Int}}}} = [[] for _ in 0:tot_deg_max]
        for (i, vec) in enumerate(indexes)
            for (idx, val) in vec
                push!(tmp[i], (0, idx...) => val)
                v = 1
                for k in 1:(tot_deg_max - i + 1)
                    v *= k//(k + i - 1 + dim)
                    push!(tmp[i + k], (k, idx...) => val * v)
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

function (ec::EmbeddedCubature{H,L,D,T})(fct::Function) where {H,L,D,T<:Real}
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
