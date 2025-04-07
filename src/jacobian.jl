# using ForwardDiff
using EcologicalNetworksDynamics
import ForwardDiff: jacobian

"""
    get_dBdt(m::Model)

Generate the function returning the vector of
species growth rate (dB/dt) given the vector of their biomass (B),
given a `Model` from EcologicalNetworksDynamics.

The output is aimed to be passed to the `jacobian` function of `Zygote`.
"""
function get_dBdt(m::Model)
    dudt = EcologicalNetworksDynamics.Internals.dudt
    (B) -> dudt(B, m._value)
end
export get_dBdt

function jacobian(m::Model, B::AbstractVector)
    jacobian(get_dBdt(m), B)
end
export jacobian
