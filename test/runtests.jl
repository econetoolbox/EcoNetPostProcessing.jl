using Test

using EcologicalNetworksDynamics
using EcoNetDynOutputs
using LinearAlgebra
show_degenerated = false # Silent warnings.

# Test files.
#include("./test-utils.jl")
include("./test-jacobian.jl")
include("./test-sensitivity.jl")
include("./doctests.jl")
