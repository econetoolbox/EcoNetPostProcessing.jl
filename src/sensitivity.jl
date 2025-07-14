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

"""
    resistance(m::Model, B::AbstractVector; response_of = :all, perturbation_on = :all, aggregated = false)

Compute the resistance of species or the entire community to a press disturbance (such as an increase in mortality).
Resistance computation is based on the sensitivity matrix.
To get resistance computed from simulated dynamics use [`resistance_simulation`](@ref).

### Main arguments

  - `m` specifies the dynamical model.
  - `B` is the vector of species equilibrium biomass.

### Keyword arguments

  - `response_of` the response_of which species is measured, by default set to `:all`.
  - `perturbation_on` which species is affected by the press, by default set to `:all`.
  - `agregated` if true return the resistance on the species selected by `response_of`,
    can be used typically to compute the resistance at the community level.
    By default set to false, so that resistance is evaluated at the species level.

See also [`sensitivity`](@ref).
"""
function resistance(
    m::Model,
    B::AbstractVector;
    response_of = :all,
    perturbation_on = :all,
    aggregated = false,
)
    s = sensitivity(m, B)
    S = m.richness
    # Process response_of.
    if response_of == :all
        rows = 1:S
    elseif eltype(response_of) <: Integer
        rows = response_of
    else
        @error "Argument `response_of` should be either set to :all or a (vector of) integer(s)."
    end
    # Process perturbation_on.
    if perturbation_on == :all
        cols = 1:S
    elseif eltype(perturbation_on) <: Integer
        cols = perturbation_on
    else
        @error "Argument `perturbation_on` should be either set to :all or a (vector of) integer(s)."
    end
    res = length(cols) > 1 ? vec(sum(s[rows, cols]; dims = 2)) : s[rows, cols]
    aggregated && (res = sum(res))
    -res
end
export resistance

"""
    resistance_simulation(
    m::Model,
    B::AbstractVector;
    mortality_increment::AbstractVector = fill(0.1, m.richness),
    response_of = :all,
    aggregated = false,
    normalized = true,

)

Compute the resistance of species or group of species to a mortality increase.
Resistance is computed from simulation outputs.
For analytical computation see [`resistance`](@ref).
For small `mortality_increment` values both function should be equivalent.

Resistance is defined as the change in biomass relative to the mortality increment, when `normalized` is true.
When `normalized` is false, the absolute change in biomass is returned.

When `aggregated` is set to true, the total change in biomass is computed, and is normalized
by the mean mortality increment.

Note that if `normalized` is true and the mortality increment of species is set to 0,
its resistance value will be set to `:undefined` because of division by 0.
To bypass this behaviour, set `:normalized` to false, and then use the normalization of your choice
on the output.

### Main arguments

  - `m` specifies the dynamical model.
  - `B` is the vector of species equilibrium biomass.

### Keyword arguments

  - `response_of` the response_of which species is measured, by default set to `:all`.
  - `mortality_increment` vector which species the increase in mortality rate for each species.
  - `agggregated` if true return the resistance on the species selected by `response_of`,
    can be used typically to compute the resistance at the community level.
    By default set to false, so that resistance is evaluated at the species level.
  - `normalized` is the change in biomass normalized by the mortality increment. By default set to `true`.

See also [`sensitivity`](@ref).
"""
function resistance_simulation(
    m::Model,
    B::AbstractVector;
    mortality_increment::AbstractVector = fill(0.1, m.richness),
    response_of = :all,
    aggregated = false,
    normalized = true,
)
    S = m.richness
    m_press = deepcopy(m)
    for i in 1:S
        m_press.d[i] += mortality_increment[i]
    end
    B_end = simulate(m_press, B, 1_000; show_degenerated = false).u[end]
    B_diff = B_end .- B
    # Process response_of.
    if response_of == :all
        rows = 1:S
    elseif eltype(response_of) <: Integer
        rows = response_of
    else
        @error "Argument `response_of` should be either set to :all or a (vector of) integer(s)."
    end
    B_diff = B_diff[rows]
    if !aggregated
        res = normalized ? B_diff ./ mortality_increment[rows] : B_diff
        res = convert(Vector{Any}, res)
        res[isnan.(res)] .= :undefined
    else
        res = normalized ? sum(B_diff) / mean(mortality_increment) : sum(B_diff)
    end
    res
end
export resistance_simulation
