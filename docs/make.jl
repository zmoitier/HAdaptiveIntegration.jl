using AdaptiveSimplexQuadrature
using Documenter

DocMeta.setdocmeta!(AdaptiveSimplexQuadrature, :DocTestSetup, :(using AdaptiveSimplexQuadrature); recursive=true)

makedocs(;
    modules=[AdaptiveSimplexQuadrature],
    authors="Zoïs Moitier and Luiz M. Faria",
    sitename="AdaptiveSimplexQuadrature.jl",
    format=Documenter.HTML(;
        canonical="https://zmoitier.github.io/AdaptiveSimplexQuadrature.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/zmoitier/AdaptiveSimplexQuadrature.jl",
    devbranch="main",
)
