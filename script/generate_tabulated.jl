using HAdaptiveIntegration: EmbeddedCubature, GenzMalik, GrundmannMoeller, embedded_cubature
using Optim
using Printf: Format, format

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

function triangle_gm()
    ec = embedded_cubature(BigFloat, GrundmannMoeller{2}(7, 5))
    print_embedded_cubature(ec, "Grundmann-Möller")
    return nothing
end

function square_gm()
    ec = embedded_cubature(BigFloat, GenzMalik{2}())
    print_embedded_cubature(ec, "Genz-Malik")
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
