using HCubature
using HAdaptiveIntegration: integrate, Cuboid
using HAdaptiveIntegration.Rule: CUBE_BE65, GenzMalik, embedded_cubature, orders
using CairoMakie
using StaticArrays
using LinearAlgebra

include("util.jl")

## Quadrature rule used for HAdaptiveIntegration (swap to GenzMalik to compare)
const QRULE = CUBE_BE65
# const QRULE = GenzMalik{3}()

function reference(f)
    I, E = hcubature(f, zeros(3), ones(3); rtol = 1.0e-12, maxevals = typemax(Int))
    return I, E
end

function run_convergence(fct)
    high, low = orders(QRULE) .+ 1
    counter_hai = Ref(0)
    counter_hc = Ref(0)
    fct_hai = (x) -> (counter_hai[] += 1; fct(x))
    fct_hc = (x) -> (counter_hc[] += 1; fct(x))

    dom = Cuboid((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))
    ec = embedded_cubature(QRULE)

    rtol_vec = [1 / 10^i for i in 1:8]
    Iref, _ = reference(fct)

    hai = (I = zeros(length(rtol_vec)), E = zeros(length(rtol_vec)), N = zeros(length(rtol_vec)))
    hc = (I = zeros(length(rtol_vec)), E = zeros(length(rtol_vec)), N = zeros(length(rtol_vec)))

    for (i, rtol) in enumerate(rtol_vec)
        counter_hai[] = 0
        I, E = integrate(fct_hai, dom; rule = ec, rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter_hai[]

        counter_hc[] = 0
        I, E = hcubature(fct_hc, zeros(3), ones(3); rtol = rtol)
        hc.I[i], hc.E[i], hc.N[i] = I, E, counter_hc[]
    end

    return hai, hc, Iref, high, low
end

## definitions
ϵ = 0.05
x₀ = SVector(1 / pi, 1 / pi, 1 / pi)
r₀ = 0.5

fct_point = (x) -> scaled_mollifier(norm(x - x₀), ϵ, 3)
fct_sphere = (x) -> scaled_mollifier(dot(x, x) - r₀^2, ϵ, 1)
fct_plane = (x) -> scaled_mollifier(x[1] - 1 / π, ϵ, 1)

println("Running convergence for point feature...")
hai_point, hc_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for sphere feature...")
hai_sphere, hc_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
println("Running convergence for plane feature...")
hai_plane, hc_plane, Iref_plane, _, _ = run_convergence(fct_plane)

fig = Figure(size = (1300, 520))

c_hai = Makie.wong_colors()[1]
c_hc = Makie.wong_colors()[2]

axes = Axis[]
legend_entries = Any[]
for (col, hai, hc, Iref, title) in (
        (1, hai_point, hc_point, Iref_point, "Point Feature"),
        (2, hai_sphere, hc_sphere, Iref_sphere, "Sphere Feature"),
        (3, hai_plane, hc_plane, Iref_plane, "Plane Feature"),
    )
    ylab = col == 1 ? "Relative error" : ""
    ax = push!(
        axes, Axis(
            fig[1, col];
            xlabel = "N (number of evaluations)",
            ylabel = ylab,
            xscale = log10,
            yscale = log10,
            title = title,
        )
    )[end]

    p1 = scatterlines!(ax, hai.N, abs.(hai.I .- Iref) ./ abs(Iref); color = c_hai, marker = :circle)
    p2 = scatterlines!(ax, hai.N, hai.E ./ abs(Iref); color = c_hai, marker = :rect)
    p3 = scatterlines!(ax, hc.N, abs.(hc.I .- Iref) ./ abs(Iref); color = c_hc, marker = :circle)
    p4 = scatterlines!(ax, hc.N, hc.E ./ abs(Iref); color = c_hc, marker = :rect)

    idx_ref = length(hai.N)
    N_ref = hai.N[idx_ref]
    e_ref = abs(hai.I[idx_ref] - Iref) / abs(Iref)
    est_ref = hai.E[idx_ref] / abs(Iref)
    p5 = lines!(ax, hai.N, e_ref .* (hai.N ./ N_ref) .^ (-high / 3); color = :black, linestyle = :dash)
    p6 = lines!(ax, hai.N, est_ref .* (hai.N ./ N_ref) .^ (-low / 3); color = :gray, linestyle = :dash)

    col == 1 && append!(legend_entries, [p1, p2, p3, p4, p5, p6])
end

linkaxes!(axes...)
Legend(
    fig[2, :],
    legend_entries,
    [
        "HAI (actual)", "HAI (estimated)",
        "HCubature (actual)", "HCubature (estimated)",
        L"\mathcal{O}(N^{-%$(high)/3})", L"\mathcal{O}(N^{-%$(low)/3})",
    ];
    orientation = :horizontal,
    framevisible = false,
)

save(joinpath(@__DIR__, "cvg_cube.png"), fig)
println("Saved cvg_cube.png")
