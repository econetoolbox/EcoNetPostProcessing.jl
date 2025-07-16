using Test

using EcologicalNetworksDynamics
using EcoNetPostProcessing
using LinearAlgebra
show_degenerated = false # Silent warnings.

# Test files.
include("./test-jacobian.jl")
include("./test-sensitivity.jl")
include("./test-robustness.jl")
include("./test-utils.jl")

# Testing doctests
println("\nRun doctests.\n")
include("./doctests.jl")
