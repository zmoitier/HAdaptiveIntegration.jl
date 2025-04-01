using StaticArrays
import Printf: Format, format

import HAdaptiveIntegration as hai

setprecision(BigFloat, 20; base=10)

function orbit(x::T, y::T, z::T) where {T}
    list = []
    for (a, b, c) in [(x, y, z), (y, x, z), (z, y, x), (x, z, y), (y, z, x), (z, x, y)]
        for (α, β, γ) in Base.product((1, -1), (1, -1), (1, -1))
            push!(list, [α * a, β * b, γ * c])
        end
    end
    sort!(list)
    return list
end

function orbit_000(w::String)
    return ([parse.(BigFloat, ("0", "0", "0"))], parse(BigFloat, w))
end

function orbit_x00(x::String, w::String)
    x_num = parse(BigFloat, x)
    z = parse(BigFloat, "0")

    list = NTuple{3,BigFloat}[]
    for α in (1, -1)
        a = α * x_num
        append!(list, [(a, z, z), (z, a, z), (z, z, a)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xx0(x::String, w::String)
    x_num = parse(BigFloat, x)
    z = parse(BigFloat, "0")

    list = NTuple{3,BigFloat}[]
    for (α, β) in Base.product((1, -1), (1, -1))
        a, b = α * x_num, β * x_num
        append!(list, [(a, b, z), (a, z, b), (z, a, b)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xy0(x::String, y::String, w::String)
    x_num, y_num = parse(BigFloat, x), parse(BigFloat, y)
    z = parse(BigFloat, "0")

    list = NTuple{3,BigFloat}[]
    for (α, β) in Base.product((1, -1), (1, -1))
        a, b = α * x_num, β * y_num
        append!(list, [(a, b, z), (b, a, z), (a, z, b), (b, z, a), (z, a, b), (z, b, a)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xxx(x::String, w::String)
    x_num = parse(BigFloat, x)
    z = parse(BigFloat, "0")

    list = NTuple{3,BigFloat}[]
    for (α, β, γ) in Base.product((1, -1), (1, -1), (1, -1))
        a, b, c = α * x_num, β * x_num, γ * x_num
        push!(list, (a, b, c))
    end

    return (list, parse(BigFloat, w))
end

function orbit_xxy(x::String, y::String, w::String)
    x_num, y_num = parse(BigFloat, x), parse(BigFloat, y)

    list = NTuple{3,BigFloat}[]
    for (α, β, γ) in Base.product((1, -1), (1, -1), (1, -1))
        a, b, c = α * x_num, β * x_num, γ * y_num
        append!(list, [(a, b, c), (a, c, b), (c, a, b)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xyz(x::String, y::String, z::String, w::String)
    x_num, y_num, z_num = parse(BigFloat, x), parse(BigFloat, y), parse(BigFloat, z)

    list = NTuple{3,BigFloat}[]
    for (α, β, γ) in Base.product((1, -1), (1, -1), (1, -1))
        a, b, c = α * x_num, β * y_num, γ * z_num
        append!(list, [(a, b, c), (b, a, c), (a, c, b), (b, c, a), (c, a, b), (c, b, a)])
    end

    return (list, parse(BigFloat, w))
end

function deorbit(orbits)
    nodes = Vector{NTuple{3,BigFloat}}()
    weights = Vector{BigFloat}()
    for (list, w) in orbits
        append!(nodes, list)
        append!(weights, fill(w, length(list)))
    end
    return (nodes=nodes, weights=weights)
end

cube = hai.cuboid([big"-1.0", big"-1.0", big"-1.0"], [big"1.0", big"1.0", big"1.0"])

# https://epubs.siam.org/doi/10.1137/0725016

# 65 points, order 9
ruleA = deorbit([
    orbit_000("3.627223234882982e-2"),
    orbit_x00("0.5964879651434033", "3.344004803960433e-1"),
    orbit_x00("0.9115074790731163", "1.056782249762152e-1"),
    orbit_xx0("0.8574202866331438", "1.052721389844229e-1"),
    orbit_xxx("0.5055319855426346", "2.134446785647350e-1"),
    orbit_xxx("0.9029552445284127", "2.932190346652714e-2"),
    orbit_xxy("0.5250000000000000", "0.9350000000000000", "8.824405047310198e-2"),
])

# 65 points, order 7
ruleN = deorbit([
    orbit_000("-1.567680589691669"),
    orbit_x00("0.5964879651434033", "7.463617511755153e-1"),
    orbit_x00("0.9115074790731163", "-2.514018470880359e-1"),
    orbit_xx0("0.8574202866331438", "-4.386770693227025e-2"),
    orbit_xxx("0.5055319855426346", "-3.526102622434519e-1"),
    orbit_xxx("0.9029552445284127", "-2.158018313159962e-2"),
    orbit_xxy("0.5250000000000000", "0.9350000000000000", "8.824405047310198e-2"),
])
ruleB = (nodes=ruleN[:nodes], weights=ruleA[:weights] - ruleN[:weights])

for (name, cbt, n) in [("BE_9_65", ruleA, 16), ("BE_7_65", ruleB, 16)]
    println(">> $name <<")

    local Φ = hai.map_to_reference(cube)
    local j = hai.abs_det_jac(hai.reference_orthotope(BigFloat, 3)) / hai.abs_det_jac(cube)

    fmt_node = Format("[\"%.$(n)e\", \"%.$(n)e\", \"%.$(n)e\"],")
    for x in cbt[:nodes]
        local v = Φ(SVector{3}(x))
        println(
            format(
                fmt_node,
                round(v[1]; sigdigits=n + 1, base=10),
                round(v[2]; sigdigits=n + 1, base=10),
                round(v[3]; sigdigits=n + 1, base=10),
            ),
        )
    end
    println()

    fmt_weight = Format("\"%.$(n)e\",")
    for w in cbt[:weights]
        println(format(fmt_weight, round(j * w; sigdigits=n + 1, base=10)))
    end
    println()
end
