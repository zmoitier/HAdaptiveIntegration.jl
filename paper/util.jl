using LinearAlgebra
using HAdaptiveIntegration: Triangle

function mollifier(r::Real)
    return exp(-r^2 / 2) / sqrt(2 * π)
end

function scaled_mollifier(r::Real, ϵ::Real, d::Int)
    return mollifier(r / ϵ) / (ϵ^d)
end

function plot_triangle(tri::Triangle)
    v1, v2, v3 = tri.vertices
    x = [v1[1], v2[1], v3[1], v1[1]]
    y = [v1[2], v2[2], v3[2], v1[2]]
    return x, y
end
