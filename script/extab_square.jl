using Base.Iterators: countfrom, partition
using LinearAlgebra
using Optim
using Printf: Format, format
using StaticArrays

import HAdaptiveIntegration as hai

include("./extab_utils.jl")

function integral_chebyshev_orthotope(d::Int, tdm::Int)
    if d ≤ 0
        return Vector{Vector{Pair{Tuple{},Rational{Int}}}}()
    end

    indexes = [[(k,) => iseven(k) ? 2//(1 - k * k) : 0//1] for k in 0:tdm]

    for n in 2:d
        tmp = [Vector{Pair{NTuple{n,Int},Rational{Int}}}() for _ in 0:tdm]
        for (td, idx_val) in zip(Iterators.countfrom(0), indexes)
            for (idx, val) in idx_val
                for k in 0:(tdm - td)
                    push!(
                        tmp[td + 1 + k],
                        (k, idx...) => val * (iseven(k) ? 2//(1 - k * k) : 0//1),
                    )
                end
            end
        end
        indexes = tmp
    end

    return indexes
end

function reduce_quarter(
    indexes::Vector{Vector{Pair{NTuple{D,Int64},Rational{Int64}}}}
) where {D}
    new = Vector{Vector{Pair{NTuple{D,Int64},Rational{Int64}}}}()
    for pairs in indexes
        tmp = [pair for pair in pairs if pair[2] ≠ 0]
        if !isempty(tmp)
            push!(new, tmp)
        end
    end
    return new
end

function increase_precision(ec)
    orbit_maps = ec[:orbit_maps]
    kh = ec[:order_high]
    kl = ec[:order_low]

    polynomials = reduce_quarter(integral_chebyshev_orthotope(2, kh))

    U, D, H, L = pack(ec[:nodes], ec[:weights_high], ec[:weights_low])
    @assert D == 2

    function F(u)
        v = zero(eltype(u))

        nodes = @view u[1:(2H)]
        wh = @view u[(2H + 1):(3H)]
        wl = @view u[(3H + 1):end]

        acos_x = acos.(nodes)

        for pairs in polynomials[1:(kl ÷ 2 + 1)]
            for (mi, val) in pairs
                v_lo, v_hi = zero(v), zero(v)

                for (i, acx, Φs) in zip(countfrom(1), partition(acos_x, 2), orbit_maps)
                    tmp = sum(prod(cos.(mi .* Φ(acx))) for Φ in Φs)
                    v_hi += wh[i] * tmp
                    if i ≤ L
                        v_lo += wl[i] * tmp
                    end
                end

                v += (v_hi - val)^2 + (v_lo - val)^2
            end
        end
        for pairs in polynomials[(kl ÷ 2 + 2):end]
            for (mi, val) in pairs
                v_hi = zero(v)

                for (w, acx, Φs) in zip(wh, partition(acos_x, 2), orbit_maps)
                    v_hi += w * sum(prod(cos.(mi .* Φ(acx))) for Φ in Φs)
                end

                v += (v_hi - val)^2
            end
        end

        return v
    end

    @show F(U)
    println()

    result = optimize(
        F, U, LBFGS(), Optim.Options(; x_tol=1e-40, g_tol=1e-50); autodiff=:forward
    )
    display(result)

    V = Optim.minimizer(result)
    return unpack(V, orbit_maps)
end

function main()
    square = hai.rectangle(BigFloat[-1, -1], BigFloat[1, 1])
    reference = hai.reference_domain(typeof(square))

    I = BigFloat[1 0; 0 1]
    R = BigFloat[0 -1; 1 0]
    S = BigFloat[-1 0; 0 1]
    ams_00 = [affine_map(I)]
    ams_x0 = [affine_map(R^i) for i in 0:3]
    ams_xx = ams_x0
    ams_xy = [[affine_map(R^i) for i in 0:3]..., [affine_map(R^i * S) for i in 0:3]...]

    # https://link.springer.com/article/10.1007/BF01389339
    # Cools-Haegemans with 25 points
    ch25 = assemble(
        9,
        7,
        (
            ams_00,
            ("0.00000000000000000000", "0.00000000000000000000"),
            "0.32363456790123456790",
            "0.67592092205970002525",
        ),
        (
            ams_x0,
            ("0.90617984593866399280", "0.00000000000000000000"),
            "0.13478507238752090312",
            "0.23092842785903867626",
        ),
        (
            ams_xx,
            ("0.53846931010568309104", "0.53846931010568309104"),
            "0.22908540022399111713",
            "0.43953907332966785983",
        ),
        (
            ams_xx,
            ("0.90617984593866399280", "0.90617984593866399280"),
            "0.56134348862428635955e-1",
            "0.82373073956971141166e-1",
        ),
        (
            ams_xy,
            ("0.90617984593866399280", "0.53846931010568309104"),
            "0.11340000000000000000",
            "0.39089597169698608216e-1",
        ),
        (
            ams_x0,
            ("0.53846931010568309104", "0.00000000000000000000"),
            "0.27228653255075070182",
            "",
        ),
    )
    # Cools-Haegemans with 21 points
    ch21 = assemble(
        7,
        5,
        (
            ams_00,
            ("0.00000000000000000000", "0.00000000000000000000"),
            "0.67592092205970002525",
            "0.61048736734452269380",
        ),
        (
            ams_x0,
            ("0.90617984593866399280", "0.00000000000000000000"),
            "0.23092842785903867626",
            "0.26364520521662754199",
        ),
        (
            ams_xx,
            ("0.53846931010568309104", "0.53846931010568309104"),
            "0.43953907332966785983",
            "0.47862867049936646804",
        ),
        (
            ams_xx,
            ("0.90617984593866399280", "0.90617984593866399280"),
            "0.82373073956971141166e-1",
            "0.10510428244787531652",
        ),
        (
            ams_xy,
            ("0.90617984593866399280", "0.53846931010568309104"),
            "0.39089597169698608216e-1",
            "",
        ),
    )

    for (cbt, name) in [(ch21, "CH21"), (ch25, "CH25")]
        println("##########")
        println("## $name")
        println("##########")
        println()

        orbit_maps = cbt[:orbit_maps]
        nodes, wh, wl = increase_precision(cbt)

        Ψ = hai.map_to_reference(square)
        j = hai.abs_det_jac(reference) / hai.abs_det_jac(square)

        fmt_node = Format("[\"%.36e\", \"%.36e\"],")
        fmt_weight = Format("\"%.36e\",")

        println(">> nodes <<")
        for (x, Φs) in zip(nodes, orbit_maps)
            for Φ in Φs
                y = Ψ(Φ(x))
                println(format(fmt_node, y[1], y[2]))
            end
        end
        println()

        println(">> weights high <<")
        for (w, Φs) in zip(wh, orbit_maps)
            z = j * w
            for _ in Φs
                println(format(fmt_weight, z))
            end
        end
        println()

        println(">> weights low <<")
        for (w, Φs) in zip(wl, orbit_maps)
            z = j * w
            for _ in Φs
                println(format(fmt_weight, z))
            end
        end
        println()
    end

    return nothing
end

main()
