module Doctest

using Documenter
import EcoNetDynOutputs

DocMeta.setdocmeta!(
    EcoNetDynOutputs,
    :DocTestSetup,
    :(using EcoNetDynOutputs);
    recursive = true,
)

doctest(EcoNetDynOutputs)

end
