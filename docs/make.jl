using Documenter
using DocumenterCitations
using EcologicalNetworksDynamics
using EcoNetPostProcessing

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

DocMeta.setdocmeta!(
    EcoNetPostProcessing,
    :DocTestSetup,
    :(using EcoNetPostProcessing, EcologicalNetworksDynamics);
    recursive=true,
)

makedocs(;
    authors="Alain Danet and Ismaël Lajaaiti",
    pages=[
        "Home" => "index.md",
        "Stability" => "stability.md",
        "Utilities" => "utils.md",
        "Functions" => "docstrings.md",
        "References" => "references.md",
    ],
    sitename="EcoNetPostProcessing.jl",
    repo="https://github.com/econetoolbox/EcoNetPostProcessing.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    modules=[EcoNetPostProcessing],
    plugins=[bib],
)

deploydocs(;
    repo="github.com/econetoolbox/EcoNetPostProcessing.jl",
    devbranch="main",
    branch="gh-pages",
)

