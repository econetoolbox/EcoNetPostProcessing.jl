module Doctest

using Documenter
import EcoNetPostProcessing

DocMeta.setdocmeta!(
    EcoNetPostProcessing,
    :DocTestSetup,
    :(using EcoNetPostProcessing, EcologicalNetworksDynamics);
    recursive = true,
)

doctest(EcoNetPostProcessing)

end
