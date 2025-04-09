"""
    get_pgr(m::Model)

Generate the function returning the vector of
species per capita growth rate (1/B x dB/dt)
given the vector of their biomass (B),
given a `Model` from EcologicalNetworksDynamics.

See also [`sensitivity`](@ref).
"""
function get_pgr(m::Model)
    pgr = EcologicalNetworksDynamics.Internals.dudt_per_capita
    (B) -> pgr(B, m._value)
end
export get_pgr

"""
    get_interaction(m::Model, B::AbstractVector)

Compute the interaction of the model `m`.
Because interactions are density-dependent,
the vectory of species biomass `B` should be specified.
`A[i, j]` is the interaction from species j to species i,
and is formally defined as the derivative of the
per capita growth of species i with respect to variation in biomass of species j.
See [Novak et al. 2016](https://doi.org/10.1146/annurev-ecolsys-032416-010215).

See also [`sensitivity`](@ref).
"""
function get_interaction(m::Model, B::AbstractVector)
    A = jacobian(get_pgr(m), B)
    Float64.(A)
end
export get_interaction

"""
    sensitivity(m::Model, B::AbstractVector)

Compute the sensitivity matrix of the model `m`.
Because interactions are density-dependent,
the vectory of species biomass `B` should be specified.
The sensitivity matrix is defined as the inverse of the interaction matrix.
See [Novak et al. 2016](https://doi.org/10.1146/annurev-ecolsys-032416-010215).
"""
function sensitivity(m::Model, B::AbstractVector) end
