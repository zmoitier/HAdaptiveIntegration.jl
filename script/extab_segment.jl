using Base.Iterators: countfrom, partition
using LinearAlgebra
using Optim
using Printf: Format, format
using StaticArrays

import HAdaptiveIntegration as hai

include("./extab_utils.jl")

function increase_precision(quad)
    orbit_maps = quad[:orbit_maps]
    kh = quad[:order_high]
    kl = quad[:order_low]

    U, D, H, L = pack(quad[:nodes], quad[:weights_high], quad[:weights_low])
    @assert D == 1

    function F(u)
        v = zero(eltype(u))

        x = @view u[1:H]
        wh = @view u[(H + 1):(2H)]
        wl = @view u[(2H + 1):end]

        acos_x = acos.(x)

        for k in 0:2:(kl - 1)
            val = 2//(1 - k^2)
            v_lo, v_hi = zero(v), zero(v)

            for (i, acx, Φs) in zip(countfrom(1), acos_x, orbit_maps)
                tmp = sum(cos(k * Φ(acx)[1]) for Φ in Φs)
                v_hi += wh[i] * tmp
                if i ≤ L
                    v_lo += wl[i] * tmp
                end
            end

            v += (v_hi - val)^2 + (v_lo - val)^2
        end
        for k in (kl + 1):2:kh
            val = 2//(1 - k^2)
            v_hi = zero(v)

            for (w, acx, Φs) in zip(wh, acos_x, orbit_maps)
                v_hi += w * sum(cos(k * Φ(acx)[1]) for Φ in Φs)
            end

            v += (v_hi - val)^2
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
    segment = hai.segment(big"-1.0", big"1.0")
    reference = hai.reference_domain(typeof(segment))

    ams_0 = [affine_map([big"1.0";;])]
    ams_x = [affine_map([big"1.0";;]), affine_map([big"-1.0";;])]

    # Values from QuadGK.jl
    # Gauss-Kronrod with 7 points
    gk7 = assemble(
        11,
        5,
        (ams_0, ("0.0000000000000000",), "0.45091653865847414", "0.8888888888888885"),
        (ams_x, ("0.7745966692414834",), "0.26848808986833345", "0.5555555555555556"),
        (ams_x, ("0.4342437493468026",), "0.40139741477596220", ""),
        (ams_x, ("0.9604912687080203",), "0.10465622602646725", ""),
    )
    # Gauss-Kronrod with 15 points
    gk15 = assemble(
        23,
        13,
        (ams_0, ("0.0000000000000000",), "0.20948214108472793", "0.41795918367346907"),
        (ams_x, ("0.4058451513773972",), "0.19035057806478559", "0.38183005050511887"),
        (ams_x, ("0.7415311855993945",), "0.14065325971552592", "0.27970539148927670"),
        (ams_x, ("0.9491079123427585",), "0.06309209262997842", "0.12948496616886981"),
        (ams_x, ("0.2077849550078984",), "0.20443294007529877", ""),
        (ams_x, ("0.5860872354676911",), "0.16900472663926788", ""),
        (ams_x, ("0.8648644233597691",), "0.10479001032225017", ""),
        (ams_x, ("0.9914553711208126",), "0.02293532201052925", ""),
    )
    # Gauss-Kronrod with 31 points
    gk31 = assemble(
        47,
        29,
        (ams_0, ("0.0000000000000000",), "0.10133000701479164", "0.20257824192556112"),
        (ams_x, ("0.2011940939974345",), "0.09917359872179185", "0.19843148532711144"),
        (ams_x, ("0.3941513470775634",), "0.09312659817082533", "0.18616100001556220"),
        (ams_x, ("0.5709721726085388",), "0.08308050282313299", "0.16626920581699406"),
        (ams_x, ("0.7244177313601701",), "0.06985412131872824", "0.13957067792615430"),
        (ams_x, ("0.8482065834104272",), "0.05348152469092807", "0.10715922046717188"),
        (ams_x, ("0.9372733924007060",), "0.03534636079137582", "0.07036604748810820"),
        (ams_x, ("0.9879925180204855",), "0.01500794732931611", "0.03075324199611720"),
        (ams_x, ("0.1011420669187175",), "0.10076984552387552", ""),
        (ams_x, ("0.2991800071531688",), "0.09664272698362363", ""),
        (ams_x, ("0.4850818636402397",), "0.08856444305621179", ""),
        (ams_x, ("0.6509967412974169",), "0.07684968075772045", ""),
        (ams_x, ("0.7904185014424660",), "0.06200956780067060", ""),
        (ams_x, ("0.8972645323440819",), "0.04458975132476482", ""),
        (ams_x, ("0.9677390756791391",), "0.02546084732671520", ""),
        (ams_x, ("0.9980022986933971",), "0.00537747987292333", ""),
    )

    for (quad, name) in [(gk7, "GK7"), (gk15, "GK15"), (gk31, "GK31")]
        println("##########")
        println("## $name")
        println("##########")
        println()

        orbit_maps = quad[:orbit_maps]
        nodes, wh, wl = increase_precision(quad)

        Ψ = hai.map_to_reference(segment)
        j = hai.abs_det_jac(reference) / hai.abs_det_jac(segment)

        fmt_node = Format("[\"%.36e\"],")
        fmt_weight = Format("\"%.36e\",")

        println(">> nodes <<")
        for (x, Φs) in zip(nodes, orbit_maps)
            for Φ in Φs
                y = SVector(Φ(x[1]))
                println(format(fmt_node, Ψ(y)[1]))
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
