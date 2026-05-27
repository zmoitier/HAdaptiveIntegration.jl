using CairoMakie
using HAdaptiveIntegration.Rule: GrundmannMoeller, embedded_cubature, orders
using HAdaptiveIntegration: Tetrahedron, integrate
using LinearAlgebra
using StaticArrays

include("util.jl")

const QRULE = GrundmannMoeller{3}(7, 5)

function reference(f)
    ec = embedded_cubature(QRULE)
    dom = Tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    I, E = integrate(f, dom; rule = ec, rtol = REFTOL, maxsubdiv = typemax(Int))
    return I, E
end

function run_convergence(fct)
    high, low = orders(QRULE) .+ 1
    counter = Ref(0)
    fct_count = (x) -> (counter[] += 1; fct(x))
    dom = Tetrahedron((0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    ec = embedded_cubature(QRULE)

    rtol_vec = [1 / 10^i for i in 1:8]
    Iref, _ = reference(fct)

    hai = (I = zeros(length(rtol_vec)), E = zeros(length(rtol_vec)), N = zeros(length(rtol_vec)))

    for (i, rtol) in enumerate(rtol_vec)
        counter[] = 0
        I, E = integrate(fct_count, dom; rule = ec, rtol = rtol)
        hai.I[i], hai.E[i], hai.N[i] = I, E, counter[]
    end

    return hai, Iref, high, low
end

fct_point, fct_sphere, fct_plane = make_features(3)

println("Running convergence for point feature...")
hai_point, Iref_point, high, low = run_convergence(fct_point)
println("Running convergence for hypersphere feature...")
hai_sphere, Iref_sphere, _, _ = run_convergence(fct_sphere)
println("Running convergence for hyperplane feature...")
hai_plane, Iref_plane, _, _ = run_convergence(fct_plane)

fig_cvg = Figure(size = (1200, 500))

c_hai = Makie.wong_colors()[1]

axes_cvg = Axis[]
legend_plots = Any[]
for (col, hai, Iref, title) in (
        (1, hai_point, Iref_point, "Point Feature"),
        (2, hai_sphere, Iref_sphere, "Hypersphere Feature"),
        (3, hai_plane, Iref_plane, "Hyperplane Feature"),
    )
    ylab = col == 1 ? "Relative error" : ""
    ax = push!(
        axes_cvg, Axis(
            fig_cvg[1, col];
            xlabel = "N (number of evaluations)",
            ylabel = ylab,
            xscale = log10,
            yscale = log10,
            title = title,
        )
    )[end]

    p1 = scatterlines!(ax, hai.N, abs.(hai.I .- Iref) ./ abs(Iref); color = c_hai, marker = :circle)
    p2 = scatterlines!(ax, hai.N, hai.E ./ abs(Iref); color = c_hai, marker = :rect, linestyle = :dash)
    idx_ref = length(hai.N)
    p3 = lines!(
        ax, hai.N, (abs(hai.I[idx_ref] - Iref) / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-high / 3);
        color = :black, linestyle = :dot, linewidth = 2,
    )
    p4 = lines!(
        ax, hai.N, (hai.E[idx_ref] / abs(Iref)) .* (hai.N ./ hai.N[idx_ref]) .^ (-low / 3);
        color = :gray, linestyle = :dot, linewidth = 2,
    )

    col == 1 && append!(legend_plots, [p1, p2, p3, p4])
end

linkaxes!(axes_cvg...)
Legend(
    fig_cvg[2, :],
    legend_plots,
    ["Actual error", "Estimated error", L"\mathcal{O}(N^{-%$(high)/3})", L"\mathcal{O}(N^{-%$(low)/3})"];
    orientation = :horizontal,
    framevisible = false,
)

save(joinpath(@__DIR__, "cvg_tetrahedron.png"), fig_cvg)
println("Saved cvg_tetrahedron.png")
