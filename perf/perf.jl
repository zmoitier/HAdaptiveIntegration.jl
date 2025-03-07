using LinearAlgebra, StaticArrays, BenchmarkTools
import HAdaptiveIntegration as HAI

function measure_perf(
    name,
    domain,
    fct;
    ec=HAI.default_embedded_cubature(domain),
    subdiv_algo=HAI.default_subdivision(domain),
)
    println(">> $name <<")

    counter = Ref(0)
    fct_count = x -> (counter[] += 1; fct(x))
    I, E = HAI.integrate(fct_count, domain, ec, subdiv_algo)
    @show I E counter[]
    println()

    bm = @benchmark HAI.integrate($fct, $domain, $ec, $subdiv_algo)
    display(bm)
    println()
    return nothing
end

function check_all()
    e = exp(1)
    fct = x -> cos(e * x[1] + prod(x))

    measure_perf("Segment", HAI.segment(0.0, 1.0), fct)
    measure_perf("Rectangle", HAI.rectangle((0.0, 0.0), (1.0, 1.0)), fct)
    measure_perf("Triangle", HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)), fct)
    measure_perf("Cuboid", HAI.cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0)), fct)
    measure_perf(
        "Tetrahedron",
        HAI.tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        fct,
    )

    return nothing
end

function triangle_subdiv(case::Int=0)
    domain = HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

    if case == 1
        x₀ = SVector(-0.01, 0.0)
        fct = x -> norm(x - x₀)^(-3)
    elseif case == 2
        fct = x -> 1 / norm(x)
    else
        e = exp(1)
        fct = x -> cos(e * sum(x) + prod(x))
    end

    measure_perf("subdivide_triangle4", domain, fct; subdiv_algo=HAI.subdivide_triangle4)
    measure_perf("subdivide_triangle2", domain, fct; subdiv_algo=HAI.subdivide_triangle2)

    return nothing
end

function triangle_rule()
    domain = HAI.triangle((0.0, 0.0), (1.0, 0.0), (0.0, 1.0))

    e = exp(1)
    fct = x -> cos(e * sum(x) + prod(x))

    measure_perf("TRIANGLE_RL19", domain, fct; ec=HAI.embedded_cubature(HAI.TRIANGLE_RL19))
    measure_perf("TRIANGLE_GM20", domain, fct; ec=HAI.embedded_cubature(HAI.TRIANGLE_GM20))

    return nothing
end
