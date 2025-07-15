using Documenter
using DocumenterCitations
using EcologicalNetworksDynamics
using EcoNetDynOutputs

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    style=:authoryear
)

DocMeta.setdocmeta!(
    EcoNetDynOutputs,
    :DocTestSetup,
    :(using EcoNetDynOutputs, EcologicalNetworksDynamics);
    recursive=true,
)

makedocs(;
    authors="Alain Danet & Ismaël Lajaaiti",
    pages=[
        "Home" => "index.md",
        "Stability" => "stability.md",
        "Utilities" => "utils.md",
        "Functions" => "docstrings.md",
        "References" => "references.md",
    ],
    sitename="EcoNetDynOutputs.jl",
    repo="https://github.com/econetoolbox/EcoNetDynOutputs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    modules=[EcoNetDynOutputs],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/econetoolbox/EcoNetDynOutputs.jl",
    devbranch="main",
    branch="gh-pages",
)

