using Documenter
using HAdaptiveIntegration
using HAdaptiveIntegration.Domain
using HAdaptiveIntegration.Rule

DocMeta.setdocmeta!(
    HAdaptiveIntegration, :DocTestSetup, :(using HAdaptiveIntegration); recursive=true
)

makedocs(;
    modules=[HAdaptiveIntegration],
    authors="Zois Moitier,Luiz M. Faria",
    repo="https://github.com/zmoitier/HAdaptiveIntegration.jl/blob/{commit}{path}#{line}",
    sitename="HAdaptiveIntegration.jl",
    format=Documenter.HTML(;
        canonical="https://zmoitier.github.io/HAdaptiveIntegration.jl"
    ),
    pages=[
        "Home" => "index.md",
        "Advanced usage" => "advanced.md",
        # "Examples and benchmarks" => "examples.md",
        "Extensions" => "extensions.md",
        "docstrings.md",
    ],
    pagesonly=true,
)

deploydocs(; repo="github.com/zmoitier/HAdaptiveIntegration.jl", devbranch="main")
