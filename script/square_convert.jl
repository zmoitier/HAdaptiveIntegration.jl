using StaticArrays
import Printf: Format, format, @printf

import HAdaptiveIntegration as hai

setprecision(BigFloat, 40; base=10)

function orbit_00(w::String)
    return ([parse.(BigFloat, ("0", "0"))], parse(BigFloat, w))
end

function orbit_x0(x::String, w::String)
    x_num = parse(BigFloat, x)
    z = parse(BigFloat, "0")

    list = NTuple{2,BigFloat}[]
    for α in (1, -1)
        a = α * x_num
        append!(list, [(a, z), (z, a)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xx(x::String, w::String)
    x_num = parse(BigFloat, x)

    list = NTuple{2,BigFloat}[]
    for (α, β) in Base.product((1, -1), (1, -1))
        a, b = α * x_num, β * x_num
        append!(list, [(a, b)])
    end

    return (list, parse(BigFloat, w))
end

function orbit_xy(x::String, y::String, w::String)
    x_num, y_num = parse(BigFloat, x), parse(BigFloat, y)

    list = NTuple{2,BigFloat}[]
    for (α, β) in Base.product((1, -1), (1, -1))
        a, b = α * x_num, β * y_num
        append!(list, [(a, b), (b, a)])
    end

    return (list, parse(BigFloat, w))
end

function deorbit(orbits)
    nodes = Vector{NTuple{2,BigFloat}}()
    weights = Vector{BigFloat}()
    for (list, w) in orbits
        append!(nodes, list)
        append!(weights, fill(w, length(list)))
    end
    return (nodes=nodes, weights=weights)
end

square = hai.rectangle(BigFloat, (-1, -1), (1, 1))

# https://link.springer.com/article/10.1007/BF01389339

# 25 points Gauss rule of order 9
ch25 = deorbit([
    orbit_00("0.32363456790123456790"),
    orbit_x0("0.90617984593866399280", "0.13478507238752090312"),
    orbit_xx("0.53846931010568309104", "0.22908540022399111713"),
    orbit_xx("0.90617984593866399280", "0.56134348862428635955e-1"),
    orbit_xy("0.90617984593866399280", "0.53846931010568309104", "0.1134"),
    orbit_x0("0.53846931010568309104", "0.27228653255075070182"),
])

# 21 points Cools Haegemans rule of order 7
ch21 = deorbit([
    orbit_00("0.67592092205970002525"),
    orbit_x0("0.90617984593866399280", "0.23092842785903867626"),
    orbit_xx("0.53846931010568309104", "0.43953907332966785983"),
    orbit_xx("0.90617984593866399280", "0.82373073956971141166e-1"),
    orbit_xy(
        "0.90617984593866399280", "0.53846931010568309104", "0.39089597169698608216e-1"
    ),
])

# 13 points Cools Haegemans rule of order 5
ch13 = deorbit([
    orbit_00("0.61048736734452269380"),
    orbit_x0("0.90617984593866399280", "0.26364520521662754199"),
    orbit_xx("0.53846931010568309104", "0.47862867049936646804"),
    orbit_xx("0.90617984593866399280", "0.10510428244787531652"),
])

n = 19
for (name, cbt) in [("CH25", ch25), ("CH21", ch21), ("CH13", ch13)]
    println(">> $name <<")

    local Φ = hai.map_to_reference(square)
    local j =
        hai.abs_det_jac(hai.reference_orthotope(BigFloat, 2)) / hai.abs_det_jac(square)

    fmt_node = Format("[\"%.$(n)e\", \"%.$(n)e\"],")
    for x in cbt[:nodes]
        local v = Φ(SVector{2}(x))
        println(
            format(
                fmt_node,
                round(v[1]; sigdigits=n + 1, base=10),
                round(v[2]; sigdigits=n + 1, base=10),
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
