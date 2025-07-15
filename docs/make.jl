using Documenter
using EcologicalNetworksDynamics
using EcoNetDynOutputs

DocMeta.setdocmeta!(
    EcoNetDynOutputs,
    :DocTestSetup,
    :(using EcoNetDynOutputs, EcologicalNetworksDynamics);
    recursive=true,
)

makedocs(;
    pages=[
        "Home" => "index.md",
        "Stability" => "stability.md",
        "Utilities" => "utils.md",
        "Functions" => "docstrings.md"
    ],
    sitename="EcoNetDynOutputs.jl",
    repo="https://github.com/econetoolbox/EcoNetDynOutputs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        assets=String[],
    ),
    modules=[EcoNetDynOutputs],
)

deploydocs(;
    repo="github.com/econetoolbox/EcoNetDynOutputs.jl",
    devbranch="main",
    branch="gh-pages",
)

