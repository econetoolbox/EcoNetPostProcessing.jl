module Doctest

using Documenter
import EcoNetDynOutputs

DocMeta.setdocmeta!(
    EcoNetDynOutputs,
    :DocTestSetup,
    :(using EcoNetDynOutputs, EcologicalNetworksDynamics);
    recursive = true,
)

doctest(EcoNetDynOutputs)

end
