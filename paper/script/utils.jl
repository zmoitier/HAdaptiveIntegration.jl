using LinearAlgebra
using StaticArrays

## Relative tolerance for reference solutions
const REFTOL = 1.0e-12

## Feature parameters shared across all example scripts
const FEATURE_ϵ = 0.05
const FEATURE_r₀ = 2 / 3
const FEATURE_c₀ = 1 / π

"""
    make_features(d) -> (fct_point, fct_sphere, fct_plane)

Return the three standard test integrands in `d` dimensions with the shared parameters
`FEATURE_ϵ`, `FEATURE_r₀`, and `FEATURE_c₀`.
"""
function make_features(d::Int)
    x₀ = SVector(ntuple(_ -> 1 / π, d)...)
    fct_point = (x) -> scaled_mollifier(norm(x - x₀), FEATURE_ϵ, d)

    y = setindex(zeros(SVector{d, Float64}), 1.0, 1)
    function fct_sphere(x)
        return scaled_mollifier(sum(abs2, x - y) - FEATURE_r₀^2, FEATURE_ϵ, 1)
    end

    fct_plane = (x) -> scaled_mollifier(x[1] - FEATURE_c₀, FEATURE_ϵ, 1)

    return fct_point, fct_sphere, fct_plane
end

function mollifier(r::Real)
    return exp(-r^2 / 2) / sqrt(2 * π)
end

function scaled_mollifier(r::Real, ϵ::Real, d::Int)
    return mollifier(r / ϵ) / (ϵ^d)
end

"""
    figure_sizes(width_fraction, height_fraction; textwidth_bp=385.89, font_size_pt=10)
    -> (width_px, height_px, font_size_px)

Return Makie figure dimensions (in px) for a figure authored to span
`width_fraction` × `textwidth_bp` by `height_fraction` × `textwidth_bp`, along
with the font size (in px).

Keyword arguments:
- `textwidth_bp`: TeX `\\textwidth` in bp (PDF/PostScript points).  Default
  `385.89` bp = the JOSS/inara template (10pt, a4paper; geometry left=1cm,
  right=1.5cm).
- `font_size_pt`: TeX font size in pt.  Default `10` pt = JOSS body font.
"""
function figure_sizes(
        width_fraction::Real,
        height_fraction::Real;
        textwidth_bp::Real = 385.89,
        font_size_pt::Real = 10
    )

    # Makie PDF export: 1 px = 1/96 in (0.75 bp); 1 pt (TeX) = 1/72.27 in.
    textwidth_px = textwidth_bp * 96 / 72      # bp → px
    font_size_px = font_size_pt * 96 / 72.27   # pt → px

    width_px = textwidth_px * width_fraction
    height_px = textwidth_px * height_fraction

    return width_px, height_px, font_size_px
end

function paper_theme(font_size::Real)
    return Theme(
        fontsize = font_size,
        figure_padding = 4,
        Axis = (
            titlesize = font_size - 2,
            xlabelsize = font_size,
            xticklabelsize = font_size - 2,
            yticklabelsize = font_size - 2,
            spinewidth = 0.6,
            xgridwidth = 0.4,
            ygridwidth = 0.4,
        ),
        Lines = (linewidth = 1.0,),
        ScatterLines = (linewidth = 1.0, markersize = 4),
        Legend = (framevisible = false, labelsize = font_size, patchsize = (16.0f0, 8.0f0)),
    )
end
