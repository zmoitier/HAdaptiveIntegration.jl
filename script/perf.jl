using LinearAlgebra, StaticArrays, BenchmarkTools
import HAdaptiveIntegration as hai

function get_fct(dim::Int, case::Int=0)
    if case == 1
        print("-- Nearly-singular --\n\n")
        x₀ = SVector{dim,Float64}([-0.1, [0.0 for _ in 2:dim]...])
        return x -> norm(x - x₀)^(-3)

    elseif case == 2
        print("-- Singular --\n\n")
        return x -> 1 / norm(x)

    else
        print("-- Regular --\n\n")
        e = exp(1)
        return x -> cos(e * x[1] + prod(x))
    end
end

function measure_perf(
    name,
    domain,
    fct;
    ec=hai.default_embedded_cubature(domain),
    subdiv_algo=hai.default_subdivision(domain),
    buffer=nothing,
)
    println(">> $name <<")

    rtol = 1e-8

    counter = Ref(0)
    fct_count = x -> (counter[] += 1; fct(x))
    I, E = hai.integrate(
        fct_count,
        domain;
        embedded_cubature=ec,
        subdiv_algo=subdiv_algo,
        buffer=buffer,
        rtol=rtol,
    )
    @show I E counter[]
    println()

    bm = @benchmark hai.integrate(
        $fct,
        $domain,
        embedded_cubature=($ec),
        subdiv_algo=($subdiv_algo),
        buffer=($buffer),
        rtol=($rtol),
    )
    display(bm)
    println()

    return I, E, counter[]
end

function check_all()
    e = exp(1)
    fct = x -> cos(e * x[1] + prod(x))

    measure_perf("Segment", hai.segment(0, 1), fct)
    measure_perf("Rectangle", hai.rectangle((0, 0), (1, 1)), fct)
    measure_perf("Triangle", hai.triangle((0, 0), (1, 0), (0, 1)), fct)
    measure_perf("Cuboid", hai.cuboid((0, 0, 0), (1, 1, 1)), fct)
    measure_perf(
        "Tetrahedron", hai.tetrahedron((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)), fct
    )

    return nothing
end

function triangle_subdiv(case::Int=0)
    domain = hai.triangle((0, 0), (1, 0), (0, 1))
    fct = get_fct(2, case)

    function subdivide_triangle2(t::hai.Triangle{T}) where {T}
        a, b, c = t.vertices
        bc = (b + c) / 2
        return (hai.Triangle{T}(bc, a, b), hai.Triangle{T}(bc, c, a))
    end

    measure_perf("subdivide_triangle", domain, fct; subdiv_algo=hai.subdivide_triangle)
    measure_perf("subdivide_triangle2", domain, fct; subdiv_algo=hai.subdivide_triangle2)

    return nothing
end

function triangle_rule(case::Int=0)
    domain = hai.triangle((0, 0), (1, 0), (0, 1))
    fct = get_fct(2, case)

    measure_perf("TRIANGLE_RL19", domain, fct; ec=hai.embedded_cubature(hai.TRIANGLE_RL19))
    measure_perf("TRIANGLE_GM19", domain, fct; ec=hai.embedded_cubature(hai.TRIANGLE_GM19))
    measure_perf(
        "TRIANGLE_GM_9_7",
        domain,
        fct;
        ec=hai.embedded_cubature(hai.GrundmannMoeller{2}(9, 7)),
    )

    return nothing
end

function triangle_duffy(case::Int=0)
    if case == 1
        println("-- Nearly-singular --")
        x₀ = SVector(-0.1, 0)
        f_tr = x -> norm(x - x₀)^(-3)
        f_sq = u -> f_tr(SVector(u[1], u[1] * u[2])) * u[1]
        I_ref = 3.51460795210773207
    elseif case == 2
        println("-- Singular --")
        f_tr = x -> cos(sum(x) + prod(x)) / norm(x)
        f_sq = u -> u[1] * f_tr(SVector(u[1], u[1] * u[2]))
        I_ref = 4.70961690636928440e-1
    else
        println("-- Regular --")
        e = exp(1)
        f_tr = x -> cos(e * sum(x) + prod(x))
        f_sq = u -> f_tr(SVector(u[1], u[1] * u[2])) * u[1]
        I_ref = -1.8612872839611775e-1
    end
    # f_sq = u -> f_tr(SVector(u[1], u[1] * u[2])) * u[1]

    # I, E = hai.integrate(
    #     f_tr,
    #     hai.triangle(BigFloat, (0, 0), (1, 0), (1, 1));
    #     embedded_cubature=hai.embedded_cubature(BigFloat, hai.GrundmannMoeller(2, 7)),
    #     maxsubdiv=2e4,
    # )
    # @show I E

    triangle = hai.triangle((0, 0), (1, 0), (1, 1))
    buffer_tr = hai.allocate_buffer(f_tr, triangle)

    square = hai.rectangle((0, 0), (1, 1))
    buffer_sq = hai.allocate_buffer(f_tr, square)

    It, Et, ct = measure_perf("Triangle", triangle, f_tr; buffer=buffer_tr)
    Is, Es, cs = measure_perf("Square", square, f_sq; buffer=buffer_sq)

    println("triangle eval-count = $ct")
    println("triangle est-err = $Et")
    println("triangle rel-err = $(abs(It/I_ref-1))")
    println("triangle abs-err = $(abs(It-I_ref))")
    println()
    println("square eval-count = $cs")
    println("square est-err = $Es")
    println("square rel-err = $(abs(Is/I_ref-1))")
    println("square abs-err = $(abs(Is-I_ref))")

    return nothing
end

function square_cut(case::Int=0)
    if case == 1
        println("-- Nearly-singular --")
        x₀ = SVector(-0.1, 0)
        fct = x -> norm(x - x₀)^(-3)
        I_ref = 8.69841482591011649
    elseif case == 2
        println("-- Singular --")
        fct = x -> 1 / norm(x)
        I_ref = 1.76274717403908604
    else
        println("-- Regular --")
        e = exp(1)
        fct = x -> cos(e * sum(x) + prod(x))
        I_ref = -3.72257456792235497e-1
    end

    # t1 = hai.triangle(BigFloat, (0, 0), (1, 0), (0, 1))
    # t2 = hai.triangle(BigFloat, (1, 1), (0, 1), (1, 0))
    # ec = hai.embedded_cubature(BigFloat, hai.GrundmannMoeller(2, 7))
    # buffer = hai.allocate_buffer(fct, t1, ec)
    # n = 2e4

    # I1, E1 = hai.integrate(fct, t1; embedded_cubature=ec, buffer=buffer, maxsubdiv=n)
    # I2, E2 = hai.integrate(fct, t2; embedded_cubature=ec, buffer=buffer, maxsubdiv=n)
    # @show I1 + I2 max(E1, E2)

    t1 = hai.triangle((0, 0), (1, 0), (0, 1))
    t2 = hai.triangle((1, 1), (0, 1), (1, 0))
    buffer_tr = hai.allocate_buffer(fct, t1)

    sq = hai.rectangle((0, 0), (1, 1))
    buffer_sq = hai.allocate_buffer(fct, sq)

    I1, E1, c1 = measure_perf("Triangle 1", t1, fct; buffer=buffer_tr)
    I2, E2, c2 = measure_perf("Triangle 2", t2, fct; buffer=buffer_tr)
    Is, Es, cs = measure_perf("Square", sq, fct; buffer=buffer_sq)

    println("triangle eval-count = $(c1+c2)")
    println("triangle est-err = $(max(E1,E2))")
    println("triangle rel-err = $(abs((I1+I2)/I_ref-1))")
    println("triangle abs-err = $(abs(I1+I2-I_ref))")
    println()
    println("square eval-count = $cs")
    println("square est-err = $Es")
    println("square rel-err = $(abs(Is/I_ref-1))")
    println("square abs-err = $(abs(Is-I_ref))")

    return nothing
end
