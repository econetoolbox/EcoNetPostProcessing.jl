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
sensitivity(m::Model, B::AbstractVector) = -inv(get_interaction(m, B))
export sensitivity

"""
    keystoneness(m::Model, B::AbstractVector)

Compute species keystoneness from a given model `m`
at the point given by the vector of species biomass `B`.
Most of the time `B` should be the vector of species equilibrium biomass.
Keystoneness of species i is defined as the sum `abs(S[i, j])` for
j different from i.
It quantifies how a change in the growth rate of species i
impacts all other species.

For a formal definition see [Li et al. 2025](https://doi.org/10.1111/ele.70086).

See also [`sensitivity`](@ref).
"""
function keystoneness(m::Model, B::AbstractVector)
    S = sensitivity(m, B)
    S_nodiag = S - Diagonal(S)
    vec(sum(abs, S_nodiag; dims = 1))
end
export keystoneness
