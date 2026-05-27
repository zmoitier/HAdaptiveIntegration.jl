using Quadmath
using HAdaptiveIntegration: integrate, Tetrahedron
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature, orders
using CairoMakie
using StaticArrays
using LinearAlgebra

include("util.jl")

function reference(f)
    T = Float64
    f_ext = (x) -> f(T.(x))
    rule = embedded_cubature(GrundmannMoeller{3}(7, 5), T)
    dom = Tetrahedron{T}((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    I, E = integrate(f_ext, dom; rtol = T(1.0e-10), rule = rule, maxsubdiv = 10^6)
    return I, E
end

function run_convergence(fct)
    high, low = orders(GrundmannMoeller{3}(7, 5)) .+ 1
    counter = Ref(0)
    fct_count = (x) -> (counter[] += 1; fct(x))
    dom = Tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))

    rtol_vec = [1 / 10^i for i in 1:8]
    Iref, Eref = reference(fct)

    hai = (
        I = zeros(length(rtol_vec)),
        E = zeros(length(rtol_vec)),
        N = zeros(length(rtol_vec)),
    )

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, dom; rtol = rtol)
        hai.I[i] = I
        hai.E[i] = E
        hai.N[i] = counter[]
    end

    return hai, Iref, high, low
end

## definitions
ϵ = 0.05
x₀ = SVector(1 / 3, 1 / 5, 1 / 7)
r₀ = 0.5

# Standardized integrands for mesh plot and convergence
fct_point = (x) -> scaled_mollifier(norm(x - x₀), ϵ, 3)
fct_curve = (x) -> scaled_mollifier(abs(norm(x) - r₀), ϵ, 3)

# Generate Convergence Plots
println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for curve feature...")
hai_curve, Iref_curve, _, _ = run_convergence(fct_curve)

fig_cvg = Figure(size = (800, 400))

# Point Subplot
ax_cvg1 = Axis(
    fig_cvg[1, 1], xlabel = "N (number of evaluations)", ylabel = "Relative error",
    xscale = log10, yscale = log10, title = "Point Feature"
)
scatterlines!(ax_cvg1, hai_point.N, abs.(hai_point.I .- Iref_point) ./ abs(Iref_point); marker = :circle, label = "Actual error")
scatterlines!(ax_cvg1, hai_point.N, hai_point.E ./ abs(Iref_point); marker = :rect, label = "Estimated error")

# dashed lines for slopes
lines!(
    ax_cvg1, hai_point.N, (abs(hai_point.I[end] - Iref_point) / abs(Iref_point)) .* (hai_point.N ./ hai_point.N[end]) .^ (-high / 3);
    color = :black, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(high)/3})"
)
lines!(
    ax_cvg1, hai_point.N, (hai_point.E[end] / abs(Iref_point)) .* (hai_point.N ./ hai_point.N[end]) .^ (-low / 3);
    color = :gray, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(low)/3})"
)

axislegend(ax_cvg1; position = :lb)

# Curve Subplot
ax_cvg2 = Axis(
    fig_cvg[1, 2], xlabel = "N (number of evaluations)",
    xscale = log10, yscale = log10, title = "Surface Feature"
)
scatterlines!(ax_cvg2, hai_curve.N, abs.(hai_curve.I .- Iref_curve) ./ abs(Iref_curve); marker = :circle, label = "Actual error")
scatterlines!(ax_cvg2, hai_curve.N, hai_curve.E ./ abs(Iref_curve); marker = :rect, label = "Estimated error")

lines!(
    ax_cvg2, hai_curve.N, (abs(hai_curve.I[end] - Iref_curve) / abs(Iref_curve)) .* (hai_curve.N ./ hai_curve.N[end]) .^ (-high / 3);
    color = :black, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(high)/2})"
)
lines!(
    ax_cvg2, hai_curve.N, (hai_curve.E[end] / abs(Iref_curve)) .* (hai_curve.N ./ hai_curve.N[end]) .^ (-low / 3);
    color = :gray, linestyle = :dash, label = L"\mathcal{O}(N^{-%$(low)/2})"
)

axislegend(ax_cvg2; position = :lb)

save("cvg_tetrahedron.png", fig_cvg)
println("Saved cvg_tetrahedron.png")
