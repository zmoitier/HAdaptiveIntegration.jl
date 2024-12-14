using HAdaptiveIntegration
using Documenter

DocMeta.setdocmeta!(
    HAdaptiveIntegration, :DocTestSetup, :(using HAdaptiveIntegration); recursive=true
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules=[HAdaptiveIntegration],
    authors="Zois Moitier,Luiz M. Faria",
    repo="https://github.com/zmoitier/HAdaptiveIntegration.jl/blob/{commit}{path}#{line}",
    sitename="HAdaptiveIntegration.jl",
    format=Documenter.HTML(;
        canonical="https://zmoitier.github.io/HAdaptiveIntegration.jl"
    ),
    pages=["Home" => "index.md"; numbered_pages],
)

deploydocs(; repo="github.com/zmoitier/HAdaptiveIntegration.jl", devbranch="main")
