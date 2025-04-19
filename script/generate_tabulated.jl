using Base.Iterators: countfrom, partition, product
using HAdaptiveIntegration:
    AbstractDomain,
    EmbeddedCubature,
    GenzMalik,
    GrundmannMoeller,
    RadonLaurie,
    abs_det_jac,
    cuboid,
    embedded_cubature,
    map_to_reference,
    rectangle,
    reference_domain,
    segment
using LinearAlgebra
using Optim
using Printf: Format, format
using StaticArrays

struct AffineMap{D,T,N}
    A::SMatrix{D,D,T,N}
    b::SVector{D,T}
end

function affine_map(A::AbstractMatrix, b::AbstractVector)
    n, m = size(A)
    @assert m == length(b)

    T = promote_type(eltype(A), eltype(b))
    return AffineMap(SMatrix{n,m,T}(A), SVector{m,T}(b))
end

function affine_map(A::AbstractMatrix)
    n, m = size(A)
    return AffineMap(SMatrix{n,m}(A), zero(SVector{m,eltype(A)}))
end

function (Φ::AffineMap)(x)
    return Φ.A * x + Φ.b
end

function (Φ::AffineMap{1,T,1})(x) where {T}
    return Φ.A[1] * first(x) + Φ.b[1]
end

function matrices_to_orbit(left::Vector{<:AbstractMatrix}, right::Vector{<:AbstractMatrix})
    return reduce(vcat, map(affine_map, A * B for (A, B) in product(left, right)))
end

function orbits_segment()
    return Dict(
        "0" => [affine_map([big"1.0";;])],
        "x" => [affine_map([big"1.0";;]), affine_map([big"-1.0";;])],
    )
end

function orbit_square()
    Id = Matrix{BigFloat}(I, 2, 2)
    T12 = BigFloat[0 1; 1 0]
    N1 = Diagonal(BigFloat[-1, 1])
    N2 = Diagonal(BigFloat[1, -1])
    N12 = Diagonal(BigFloat[-1, -1])

    return Dict(
        "00" => [affine_map(Id)],
        "x0" => matrices_to_orbit([Id, T12], [Id, N1]),
        "xx" => affine_map.([Id, N1, N2, N12]),
        "xy" => matrices_to_orbit([Id, T12], [Id, N1, N2, N12]),
    )
end

function orbit_cube()
    Id = Id = Matrix{BigFloat}(I, 3, 3)
    T12 = BigFloat[0 1 0; 1 0 0; 0 0 1]
    T13 = BigFloat[0 0 1; 0 1 0; 1 0 0]
    T23 = BigFloat[1 0 0; 0 0 1; 0 1 0]
    C123 = BigFloat[0 1 0; 0 0 1; 1 0 0]
    C132 = BigFloat[0 0 1; 1 0 0; 0 1 0]
    N1 = Diagonal(BigFloat[-1, 1, 1])
    N2 = Diagonal(BigFloat[1, -1, 1])
    N3 = Diagonal(BigFloat[1, 1, -1])
    N12 = Diagonal(BigFloat[-1, -1, 1])
    N13 = Diagonal(BigFloat[-1, 1, -1])
    N23 = Diagonal(BigFloat[1, -1, -1])
    N123 = Diagonal(BigFloat[-1, -1, -1])

    orbits = Dict("000" => affine_map.([Id]))
    orbits["x00"] = matrices_to_orbit([Id, T12, T13], [Id, N1])
    orbits["xx0"] = matrices_to_orbit([Id, T23, T13], [Id, N1, N2, N12])
    orbits["xxx"] = affine_map.([Id, N1, N2, N3, N12, N13, N23, N123])
    orbits["xy0"] = matrices_to_orbit([Id, T12, T23, C132, C123, T13], [Id, N1, N2, N12])
    orbits["xxy"] = matrices_to_orbit([Id, T23, T13], [Id, N1, N2, N3, N12, N13, N23, N123])
    orbits["xyz"] = matrices_to_orbit(
        [Id, T12, T23, C132, C123, T13], [Id, N1, N2, N3, N12, N13, N23, N123]
    )

    return orbits
end

function assemble(
    order_high::Int,
    order_low::Int,
    points::Tuple{Vector{AffineMap{D,T,N}},NTuple{D,String},String,String}...;
    subtraction::Bool=false,
) where {D,T,N}
    orbits = Vector{Vector{AffineMap{D,T,N}}}()
    nodes = Vector{SVector{D,T}}()
    weights_high = Vector{T}()
    weights_low = Vector{T}()

    for (maps, point, wh, wl) in points
        push!(orbits, maps)
        push!(nodes, SVector{D}(parse.(T, point)))
        push!(weights_high, parse(T, wh))
        if wl != ""
            w = subtraction ? weights_high[end] - parse(T, wl) : parse(T, wl)
            push!(weights_low, w)
        end
    end

    return (
        orbits=orbits,
        nodes=nodes,
        weights_high=weights_high,
        weights_low=weights_low,
        order_high=order_high,
        order_low=order_low,
    )
end

# [ nodes[1]..., ..., nodes[end]..., weights_high..., weights_low... ]
function pack(
    nodes::Vector{SVector{D,T}}, weights_high::Vector{T}, weights_low::Vector{T}
) where {D,T}
    U = Vector{T}()
    for node in nodes
        append!(U, node)
    end
    append!(U, weights_high)
    append!(U, weights_low)
    return U, D, length(weights_high), length(weights_low)
end

function unpack(U::Vector{T}, orbits::Vector{Vector{AffineMap{D,T,N}}}) where {D,T,N}
    H = length(orbits)
    L = length(U) - H * (D + 1)

    nodes = Vector{SVector{D,T}}()
    weights_high = Vector{T}()
    weights_low = Vector{T}()
    for (i, orbit) in enumerate(orbits)
        node = U[((i - 1) * D + 1):(i * D)]

        for Φ in orbit
            push!(nodes, SVector{D}(Φ(node)))
            push!(weights_high, U[H * D + i])
            if i ≤ L
                push!(weights_low, U[H * D + H + i])
            end
        end
    end

    return EmbeddedCubature(nodes, weights_high, weights_low)
end

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

function remove_odd(
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

function increase_precision_orthotope(ec, ::Val{D}) where {D}
    orbits = ec[:orbits]
    kh = ec[:order_high]
    kl = ec[:order_low]

    polynomials = remove_odd(integral_chebyshev_orthotope(D, kh))

    U, dim, H, L = pack(ec[:nodes], ec[:weights_high], ec[:weights_low])
    @assert D == dim

    function F(u)
        v = zero(eltype(u))

        nodes = @view u[1:(D * H)]
        wh = @view u[(D * H + 1):((D + 1) * H)]
        wl = @view u[((D + 1) * H + 1):end]

        acos_x = acos.(nodes)

        for pairs in polynomials[1:(kl ÷ 2 + 1)]
            for (mi, ref) in pairs
                vₗ, vₕ = zero(v), zero(v)

                for (i, acx, Φs) in zip(countfrom(1), partition(acos_x, D), orbits)
                    val_orbit = sum(prod(cos.(mi .* Φ(acx))) for Φ in Φs)
                    vₕ += wh[i] * val_orbit
                    if i ≤ L
                        vₗ += wl[i] * val_orbit
                    end
                end

                v += (vₕ - ref)^2 + (vₗ - ref)^2
            end
        end
        for pairs in polynomials[(kl ÷ 2 + 2):end]
            for (mi, ref) in pairs
                vₕ = zero(v)

                for (w, acx, Φs) in zip(wh, partition(acos_x, D), orbits)
                    vₕ += w * sum(prod(cos.(mi .* Φ(acx))) for Φ in Φs)
                end

                v += (vₕ - ref)^2
            end
        end

        return v
    end

    @show F(U)
    println()

    result = optimize(
        F, U, LBFGS(), Optim.Options(; x_reltol=1e-40, g_tol=1e-50); autodiff=:forward
    )
    display(result)

    V = Optim.minimizer(result)
    return unpack(V, orbits)
end

function to_zero(node::SVector{D,T}, tol::T) where {D,T}
    result = Vector{T}()
    for x in node
        push!(result, abs(x) ≤ tol ? zero(T) : x)
    end
    return SVector{D,T}(result)
end

function to_reference(ec::EmbeddedCubature{D,T}, dom::AbstractDomain{D,T}) where {D,T}
    λ = abs_det_jac(reference_domain(typeof(dom))) / abs_det_jac(dom)
    Ψ = map_to_reference(dom)

    nodes = Vector{SVector{D,T}}()
    for node in ec.nodes
        push!(nodes, Ψ(to_zero(node, T(1e-40))))
    end

    return EmbeddedCubature(nodes, λ .* ec.weights_high, λ .* ec.weights_low)
end

function print_embedded_cubature(ec::EmbeddedCubature{D,T}, name::String) where {D,T}
    fmt_node = Format("[" * join(fill("\"%.36e\"", D), ", ") * "],")
    fmt_weight = Format("\"%.36e\",")

    println("="^64)
    println("# $name")
    println()
    println("Number of nodes: ", length(ec.nodes))
    println()

    println("## nodes")
    for x in ec.nodes
        println(format(fmt_node, x...))
    end
    println()

    println("## weights_high")
    for w in ec.weights_high
        println(format(fmt_weight, w))
    end
    println()

    println("## weights_low")
    for w in ec.weights_low
        println(format(fmt_weight, w))
    end
    println()

    return nothing
end

function segment_gk7()
    orbit = orbits_segment()

    # Values from QuadGK.jl
    quad_gk = assemble(
        11,
        5,
        (orbit["0"], ("0.0000000000000000",), "0.45091653865847414", "0.8888888888888885"),
        (orbit["x"], ("0.7745966692414834",), "0.26848808986833345", "0.5555555555555556"),
        (orbit["x"], ("0.4342437493468026",), "0.40139741477596220", ""),
        (orbit["x"], ("0.9604912687080203",), "0.10465622602646725", ""),
    )
    domain = segment(big"-1.0", big"1.0")

    ec = to_reference(increase_precision_orthotope(quad_gk, Val(1)), domain)
    print_embedded_cubature(ec, "Gauss-Kronrod")
    return nothing
end

function segment_gk15()
    orbit = orbits_segment()

    # Values from QuadGK.jl
    quad_gk = assemble(
        23,
        13,
        (orbit["0"], ("0.0000000000000000",), "0.20948214108472793", "0.41795918367346907"),
        (orbit["x"], ("0.4058451513773972",), "0.19035057806478559", "0.38183005050511887"),
        (orbit["x"], ("0.7415311855993945",), "0.14065325971552592", "0.27970539148927670"),
        (orbit["x"], ("0.9491079123427585",), "0.06309209262997842", "0.12948496616886981"),
        (orbit["x"], ("0.2077849550078984",), "0.20443294007529877", ""),
        (orbit["x"], ("0.5860872354676911",), "0.16900472663926788", ""),
        (orbit["x"], ("0.8648644233597691",), "0.10479001032225017", ""),
        (orbit["x"], ("0.9914553711208126",), "0.02293532201052925", ""),
    )
    domain = segment(big"-1.0", big"1.0")

    ec = to_reference(increase_precision_orthotope(quad_gk, Val(1)), domain)
    print_embedded_cubature(ec, "Gauss-Kronrod")
    return nothing
end

function segment_gk31()
    orbit = orbits_segment()

    # Values from QuadGK.jl
    quad_gk = assemble(
        47,
        29,
        (orbit["0"], ("0.0000000000000000",), "0.10133000701479164", "0.20257824192556112"),
        (orbit["x"], ("0.2011940939974345",), "0.09917359872179185", "0.19843148532711144"),
        (orbit["x"], ("0.3941513470775634",), "0.09312659817082533", "0.18616100001556220"),
        (orbit["x"], ("0.5709721726085388",), "0.08308050282313299", "0.16626920581699406"),
        (orbit["x"], ("0.7244177313601701",), "0.06985412131872824", "0.13957067792615430"),
        (orbit["x"], ("0.8482065834104272",), "0.05348152469092807", "0.10715922046717188"),
        (orbit["x"], ("0.9372733924007060",), "0.03534636079137582", "0.07036604748810820"),
        (orbit["x"], ("0.9879925180204855",), "0.01500794732931611", "0.03075324199611720"),
        (orbit["x"], ("0.1011420669187175",), "0.10076984552387552", ""),
        (orbit["x"], ("0.2991800071531688",), "0.09664272698362363", ""),
        (orbit["x"], ("0.4850818636402397",), "0.08856444305621179", ""),
        (orbit["x"], ("0.6509967412974169",), "0.07684968075772045", ""),
        (orbit["x"], ("0.7904185014424660",), "0.06200956780067060", ""),
        (orbit["x"], ("0.8972645323440819",), "0.04458975132476482", ""),
        (orbit["x"], ("0.9677390756791391",), "0.02546084732671520", ""),
        (orbit["x"], ("0.9980022986933971",), "0.00537747987292333", ""),
    )
    domain = segment(big"-1.0", big"1.0")

    ec = to_reference(increase_precision_orthotope(quad_gk, Val(1)), domain)
    print_embedded_cubature(ec, "Gauss-Kronrod")
    return nothing
end

function triangle_gm()
    ec = embedded_cubature(BigFloat, GrundmannMoeller{2}(7, 5))
    print_embedded_cubature(ec, "Grundmann-Möller")
    return nothing
end

function triangle_rl()
    ec = embedded_cubature(BigFloat, RadonLaurie())
    print_embedded_cubature(ec, "Radon-Laurie")
    return nothing
end

function square_gm()
    ec = embedded_cubature(BigFloat, GenzMalik{2}())
    print_embedded_cubature(ec, "Genz-Malik")
    return nothing
end

function square_ch21()
    orbits = orbit_square()

    # https://link.springer.com/article/10.1007/BF01389339
    # Cools-Haegemans with 21 points
    cbt_ch = assemble(
        7,
        5,
        (
            orbits["00"],
            ("0.00000000000000000000", "0.00000000000000000000"),
            "0.67592092205970002525",
            "0.61048736734452269380",
        ),
        (
            orbits["x0"],
            ("0.90617984593866399280", "0.00000000000000000000"),
            "0.23092842785903867626",
            "0.26364520521662754199",
        ),
        (
            orbits["xx"],
            ("0.53846931010568309104", "0.53846931010568309104"),
            "0.43953907332966785983",
            "0.47862867049936646804",
        ),
        (
            orbits["xx"],
            ("0.90617984593866399280", "0.90617984593866399280"),
            "0.82373073956971141166e-1",
            "0.10510428244787531652",
        ),
        (
            orbits["xy"],
            ("0.90617984593866399280", "0.53846931010568309104"),
            "0.39089597169698608216e-1",
            "",
        ),
    )
    domain = rectangle((big"-1.0", big"-1.0"), (big"1.0", big"1.0"))

    ec = to_reference(increase_precision_orthotope(cbt_ch, Val(2)), domain)
    print_embedded_cubature(ec, "Cools-Haegemans")
    return nothing
end

function square_ch25()
    orbits = orbit_square()

    # https://link.springer.com/article/10.1007/BF01389339
    # Cools-Haegemans with 25 points
    cbt_ch = assemble(
        9,
        7,
        (
            orbits["00"],
            ("0.00000000000000000000", "0.00000000000000000000"),
            "0.32363456790123456790",
            "0.67592092205970002525",
        ),
        (
            orbits["x0"],
            ("0.90617984593866399280", "0.00000000000000000000"),
            "0.13478507238752090312",
            "0.23092842785903867626",
        ),
        (
            orbits["xx"],
            ("0.53846931010568309104", "0.53846931010568309104"),
            "0.22908540022399111713",
            "0.43953907332966785983",
        ),
        (
            orbits["xx"],
            ("0.90617984593866399280", "0.90617984593866399280"),
            "0.56134348862428635955e-1",
            "0.82373073956971141166e-1",
        ),
        (
            orbits["xy"],
            ("0.90617984593866399280", "0.53846931010568309104"),
            "0.11340000000000000000",
            "0.39089597169698608216e-1",
        ),
        (
            orbits["x0"],
            ("0.53846931010568309104", "0.00000000000000000000"),
            "0.27228653255075070182",
            "",
        ),
    )
    domain = rectangle((big"-1.0", big"-1.0"), (big"1.0", big"1.0"))

    ec = to_reference(increase_precision_orthotope(cbt_ch, Val(2)), domain)
    print_embedded_cubature(ec, "Cools-Haegemans")
    return nothing
end

function tetrahedron_gm()
    ec = embedded_cubature(BigFloat, GrundmannMoeller{3}(7, 5))
    print_embedded_cubature(ec, "Grundmann-Möller")
    return nothing
end

function cube_gm()
    ec = embedded_cubature(BigFloat, GenzMalik{3}())
    print_embedded_cubature(ec, "Genz-Malik")
    return nothing
end

function cube_be65()
    orbits = orbit_cube()

    # https://epubs.siam.org/doi/10.1137/0725016
    # Berntsen-Espelid with 65 points
    cbt_be = assemble(
        9,
        7,
        (
            orbits["000"],
            ("0.0000000000000000", "0.0000000000000000", "0.0000000000000000"),
            "3.627223234882982e-2",
            "-1.567680589691669",
        ),
        (
            orbits["x00"],
            ("0.5964879651434033", "0.0000000000000000", "0.0000000000000000"),
            "3.344004803960433e-1",
            "7.463617511755153e-1",
        ),
        (
            orbits["x00"],
            ("0.9115074790731163", "0.0000000000000000", "0.0000000000000000"),
            "1.056782249762152e-1",
            "-2.514018470880359e-1",
        ),
        (
            orbits["xx0"],
            ("0.8574202866331438", "0.8574202866331438", "0.0000000000000000"),
            "1.052721389844229e-1",
            "-4.386770693227025e-2",
        ),
        (
            orbits["xxx"],
            ("0.5055319855426346", "0.5055319855426346", "0.5055319855426346"),
            "2.134446785647350e-1",
            "-3.526102622434519e-1",
        ),
        (
            orbits["xxx"],
            ("0.9029552445284127", "0.9029552445284127", "0.9029552445284127"),
            "2.932190346652714e-2",
            "-2.158018313159962e-2",
        ),
        (
            orbits["xxy"],
            ("0.5250000000000000", "0.5250000000000000", "0.9350000000000000"),
            "8.824405047310198e-2",
            "",
        );
        subtraction=true,
    )
    domain = cuboid((big"-1.0", big"-1.0", big"-1.0"), (big"1.0", big"1.0", big"1.0"))

    ec = to_reference(increase_precision_orthotope(cbt_be, Val(3)), domain)
    print_embedded_cubature(ec, "Berntsen-Espelid")
    return nothing
end

function cube_be115()
    orbits = orbit_cube()

    # https://epubs.siam.org/doi/10.1137/0725016
    # Berntsen-Espelid with 115 points
    cbt_be = assemble(
        11,
        9,
        (
            orbits["000"],
            ("0.0000000000000000", "0.0000000000000000", "0.0000000000000000"),
            "1.324609663263692e-1",
            "7.025853877436244e-1",
        ),
        (
            orbits["x00"],
            ("0.9920271704660597", "0.0000000000000000", "0.0000000000000000"),
            "1.088626296477019e-2",
            "3.203268015578803e-2",
        ),
        (
            orbits["x00"],
            ("0.7582845050347407", "0.0000000000000000", "0.0000000000000000"),
            "2.157141671884503e-1",
            "-3.355925410762035e-2",
        ),
        (
            orbits["x00"],
            ("0.2447487465401134", "0.0000000000000000", "0.0000000000000000"),
            "4.932409085050268e-2",
            "-1.388110650417988e-1",
        ),
        (
            orbits["xx0"],
            ("0.9987344998351400", "0.9987344998351400", "0.0000000000000000"),
            "1.687648683985235e-3",
            "-1.512181927222183e-2",
        ),
        (
            orbits["xx0"],
            ("0.7793703685672423", "0.7793703685672423", "0.0000000000000000"),
            "1.346468564512807e-1",
            "1.922707982244711e-2",
        ),
        (
            orbits["xxx"],
            ("0.9999698993088767", "0.9999698993088767", "0.9999698993088767"),
            "1.750145884600386e-3",
            "-1.916070627476280e-3",
        ),
        (
            orbits["xxx"],
            ("0.7902637224771788", "0.7902637224771788", "0.7902637224771788"),
            "7.752336383837454e-2",
            "-1.854652834553552e-2",
        ),
        (
            orbits["xxx"],
            ("0.4403396687650737", "0.4403396687650737", "0.4403396687650737"),
            "2.461864902770251e-1",
            "3.727029226099722e-2",
        ),
        (
            orbits["xxy"],
            ("0.4378478459006862", "0.4378478459006862", "0.9549373822794593"),
            "6.797944868483039e-2",
            "-1.604480434502482e-2",
        ),
        (
            orbits["xxy"],
            ("0.9661093133630747", "0.9661093133630747", "0.4577105877763134"),
            "1.419962823300713e-2",
            "",
        );
        subtraction=true,
    )
    domain = cuboid((big"-1.0", big"-1.0", big"-1.0"), (big"1.0", big"1.0", big"1.0"))

    ec = to_reference(increase_precision_orthotope(cbt_be, Val(3)), domain)
    print_embedded_cubature(ec, "Berntsen-Espelid")
    return nothing
end
