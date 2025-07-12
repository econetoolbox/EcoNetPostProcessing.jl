using Documenter
using EcologicalNetworksDynamics
using EcoNetDynOutputs

makedocs(;
    sitename = "EcoNetDynOutputs",
    format = Documenter.HTML(),
    modules = [EcoNetDynOutputs],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
