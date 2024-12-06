using AdaptiveSimplexQuadrature
using Documenter

DocMeta.setdocmeta!(
    AdaptiveSimplexQuadrature,
    :DocTestSetup,
    :(using AdaptiveSimplexQuadrature);
    recursive = true,
)

const page_rename = Dict("developer.md" => "Developer docs") # Without the numbers
const numbered_pages = [
    file for file in readdir(joinpath(@__DIR__, "src")) if
    file != "index.md" && splitext(file)[2] == ".md"
]

makedocs(;
    modules = [AdaptiveSimplexQuadrature],
    authors = "Zois Moitier,Luiz M. Faria",
    repo = "https://github.com/zmoitier/AdaptiveSimplexQuadrature.jl/blob/{commit}{path}#{line}",
    sitename = "AdaptiveSimplexQuadrature.jl",
    format = Documenter.HTML(;
        canonical = "https://zmoitier.github.io/AdaptiveSimplexQuadrature.jl",
    ),
    pages = ["Home" => "index.md"; numbered_pages],
)

deploydocs(; repo = "github.com/zmoitier/AdaptiveSimplexQuadrature.jl", devbranch = "main")
