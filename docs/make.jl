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
    sitename="EcoNetDynOutputs",
    repo=Remotes.GitHub("econetoolbox", "EcoNetDynOutputs.jl"),
    format=Documenter.HTML(;
        canonical="https://github.com/econetoolbox/EcoNetDynOutputs.jl",
    ),
    modules=[EcoNetDynOutputs],
    pages=["Home" => "index.md", "Functions" => "docstrings.md"],
)

deploydocs(;
    repo="github.com/econetoolbox/EcoNetDynOutputs.jl",
    target="build", # this is where Vitepress stores its output
    devbranch="main",
    branch="gh-pages",
    push_preview=true,
)

