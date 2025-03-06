using BenchmarkTools
import HAdaptiveIntegration as HAI

function measure_perf(name, domain, fct)
    print(">> $name <<\n\n")

    ec = HAI.default_embedded_cubature(domain)
    I, E = HAI.integrate(fct, domain, ec)
    @show I E

    bm = @benchmark HAI.integrate($fct, $domain, $ec)
    display(bm)

    println()
    return nothing
end

e = exp(1)
fct = x -> cos(sum(x))
# measure_perf("Segment", HAI.segment(0.0, 1.0), fct)
# measure_perf("Rectangle", HAI.rectangle((0.0, 0.0), (1.0, 1.0)), fct)
# measure_perf("Triangle", HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)), fct)
# measure_perf("Cuboid", HAI.cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)), fct)

measure_perf(
    "Tetrahedron",
    HAI.tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    x -> cos(e * x[1] + x[2] + π * x[3]),
)

measure_perf(
    "Tetrahedron",
    HAI.tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
    x -> cos(2.718281828459045 * x[1] + x[2] + π * x[3]),
)