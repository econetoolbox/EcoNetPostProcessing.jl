# using ForwardDiff
using EcologicalNetworksDynamics
import ForwardDiff: jacobian

"""
    get_dBdt(m::Model)

Generate the function returning the vector of
species growth rate (dB/dt) given the vector of their biomass (B),
given a `Model` from EcologicalNetworksDynamics.

The output is aimed to be passed to the `jacobian` function of `ForwardDiff`.
For more information see [ForwardDiff documentation](https://juliadiff.org/ForwardDiff.jl/stable/user/api/#Jacobians-of-f(x::AbstractArray)::AbstractArray).

See also [`jacobian`](@ref).
"""
function get_dBdt(m::Model)
    dudt = EcologicalNetworksDynamics.Internals.dudt
    (B) -> dudt(B, m._value)
end
export get_dBdt

"""
    jacobian(m::Model, B::AbstractVector)

Compute the jacobian of the system specied by the model `B`.
The jacobian is evaluated in `B` which gives species biomass.
"""
function jacobian(m::Model, B::AbstractVector)
    j = jacobian(get_dBdt(m), B)
    Float64.(j)
end
export jacobian
