using LinearAlgebra
using StaticArrays
using HAdaptiveIntegration: Triangle, Rectangle

## Relative tolerance for reference solutions
const REFTOL = 1.0e-12

## Feature parameters shared across all example scripts
const FEATURE_ϵ = 0.05
const FEATURE_r₀ = 0.5
const FEATURE_c₀ = 1 / π

"""
    make_features(d) -> (fct_point, fct_sphere, fct_plane)

Return the three standard test integrands in `d` dimensions with the shared parameters
`FEATURE_ϵ`, `FEATURE_r₀`, and `FEATURE_c₀`.
"""
function make_features(d::Int)
    x₀ = SVector(ntuple(_ -> 1 / π, d)...)
    fct_point = (x) -> scaled_mollifier(norm(x - x₀), FEATURE_ϵ, d)
    fct_sphere = (x) -> scaled_mollifier(dot(x, x) - FEATURE_r₀^2, FEATURE_ϵ, 1)
    fct_plane = (x) -> scaled_mollifier(x[1] - FEATURE_c₀, FEATURE_ϵ, 1)
    return fct_point, fct_sphere, fct_plane
end

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

function plot_rectangle(rect::Rectangle)
    l, h = rect.corners[1], rect.corners[2]
    x = [l[1], h[1], h[1], l[1], l[1]]
    y = [l[2], l[2], h[2], h[2], l[2]]
    return x, y
end
