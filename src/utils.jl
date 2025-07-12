#=
Utility functions for functions in measures
=#
"""
    extract_last_timesteps(solution; idxs = nothing, quiet = false, kwargs...)

Returns the biomass matrix of species x time over the `last` timesteps.

# Arguments

  - `last`: the number of last timesteps to consider. A percentage can also be also be
    provided as a `String` ending by `%`. Defaulted to 1.
  - `idxs`: vector of species indexes or names. Set to `nothing` by default.
  - `quiet`: ignores warning issue while extracting timesteps before last species extinction

See [`richness`](@ref) for the other arguments. If `idxs` is an integer,
it returns a vector of the species biomass instead of a matrix.

# Examples

```jldoctest
julia> fw = Foodweb([0 0; 1 0])
       m = default_model(fw)
       B0, t_end = [1, 1], 1_000
       sol = simulate(m, B0, t_end);

julia> last = extract_last_timesteps(sol; last = 1, idxs = [2, 1]);
       last ≈ sol.u[end][[2, 1]]
true

julia> last = extract_last_timesteps(sol; last = 1, idxs = ["s2", "s1"]);
       last ≈ sol.u[end][[2, 1]]
true

julia> last2 = extract_last_timesteps(sol; last = 1, idxs = [2])
       last2 ≈ sol.u[end][[2]]
true

julia> last2 = extract_last_timesteps(sol; last = 1, idxs = "s2")
       last2 ≈ sol.u[end][[2]]
true
```
"""
function extract_last_timesteps(solution; idxs = nothing, quiet = false, kwargs...)
    last = process_last_timesteps(solution; quiet, kwargs...)
    out = solution[:, (end-(last-1)):end]

    # Extract species indices:
    idxs = process_idxs(solution; idxs)
    out = if isempty(out)
        zeros((length(idxs), 0))
    else
        out[idxs, :]
    end
    quiet || check_last_extinction(solution; idxs, last)

    deepcopy(out)
end
export extract_last_timesteps

"""
    process_idxs(solution; idxs = nothing)

Check and sanitize the species indices or names provided (`idxs`). Used in
[`extract_last_timesteps`](@ref) and [`living_species`](@ref).
"""
function process_idxs(solution; idxs = nothing)
    sp = get_model(solution).species.names
    # Convert from EcologicalNetworksDynamics.SpeciesNames:
    sp = string.(sp)

    if isnothing(idxs)
        idxs = sp
    end

    if idxs isa AbstractString || idxs isa Integer
        idxs = [idxs] # Convert plain values into singleton collections.
    end

    # Handle species names:
    if eltype(idxs) <: AbstractString
        idxs_in = indexin(idxs, sp)
        # Handle missing species
        absent_sp = isnothing.(idxs_in)
        any(absent_sp) && throw(
            ArgumentError("Species $(idxs[absent_sp]) are not found in the network. \
                           Any mispelling?"),
        )
        # Get the index of the species names in solution matrix
        idxs = idxs_in[.!absent_sp]
        # remove Union{nothing, ...}
        idxs = something.(idxs)
    elseif eltype(idxs) <: Integer
        check_bound_idxs = 1 .<= idxs .<= length(sp)
        absent_sp = idxs[findall(.!check_bound_idxs)]
        all(check_bound_idxs) || throw(
            ArgumentError(
                "Cannot extract idxs $(absent_sp) when there are $(length(sp)) species.",
            ),
        )
    else
        throw(ArgumentError("`idxs` should be a vector of integers (species indices) \
                             or strings (species names)"))
    end

    idxs
end
export process_idxs

function process_last_timesteps(solution; last = 1, quiet = false)

    n_timesteps = length(solution.t)

    if last isa AbstractString
        endswith(last, "%") ||
            throw(ArgumentError("The `last` argument, when given as a string, \
                                 should end with character '%'"))
        perc = parse(Float64, last[1:(end-1)])
        is_valid_perc = 0.0 < perc <= 100.0
        is_valid_perc ||
            throw(ArgumentError("Cannot extract $(perc)% of the solution's timesteps: \
                                 0% < `last` <= 100% must hold."))
        last = round(Int, n_timesteps * perc / 100)
        (last > 0 || quiet) ||
            @warn("$perc% of $n_timesteps timesteps correspond to $last output lines: \
                   an empty table has been extracted.")
    elseif last isa Integer
        last > 0 || throw(ArgumentError("Cannot extract $last timesteps. \
                                         `last` should be a positive integer."))
    elseif last isa Float64
        throw(ArgumentError("Cannot extract `last` from a floating point number. \
                             Did you mean \"$last%\"?"))
    else
        throw(ArgumentError("Cannot extract timesteps with `last=$last` \
                             of type $(typeof(last)). \
                             `last` should be a positive integer \
                             or a string representing a percentage."))
    end

    last > n_timesteps && throw(
        ArgumentError("Cannot extract $last timesteps from a trajectory solution \
                       with only $(n_timesteps) timesteps. \
                       Consider decreasing the `last` argument value \
                       and/or specifying it as a percentage instead (e.g. `\"10%\"`)."),
    )

    last
end
export process_last_timesteps

function check_last_extinction(solution; idxs = nothing, last = 1)
    ext = get_extinction_timesteps(solution; idxs)
    ext_t = ext.extinction_timestep
    n_timesteps = length(solution.t)
    check_last_extinction(n_timesteps; t = ext_t, species = ext.species, last)

end
function check_last_extinction(n_timesteps::Integer; t, species, last)
    extinct = t .!== nothing
    if any(extinct)
        check_last = findall(>(n_timesteps - (last - 1)), t[extinct])
        sp = species[extinct][check_last]
        ts = t[extinct][check_last]
        max = n_timesteps - (maximum(t[extinct]) - 1)
        isempty(check_last) ||
            @warn("With `last` = $last, a table has been extracted with the species $sp, \
                   that went extinct at timesteps = $ts. \
                   Set `last` <= $max to get rid of them.")
    end
end

function get_extinction_timesteps(solution; idxs = nothing)
    idxs = process_idxs(solution; idxs)
    sp = get_model(solution).species.names # Snippet to convert to function?
    # Convert from EcologicalNetworksDynamics.SpeciesNames:
    sp = string.(sp)[idxs]
    ext_t = findfirst.(isequal(0), eachrow(solution[idxs, :]))
    extinct = ext_t .!== nothing
    (
        species = sp[extinct],
        idxs = idxs[extinct],
        extinction_timestep = something.(ext_t[extinct]),
    )
end

function get_extinction_timesteps(m::AbstractVector; threshold = 0)
    findfirst(x -> x <= threshold, m)
end

"""
    get_alive_species(solution; idxs = nothing, threshold = 0)

Returns a tuple with species having a biomass above `threshold` at the end of a simulation.

# Examples

```jldoctest
julia> foodweb = Foodweb([0 0; 0 0]);
       m = default_model(foodweb);
       sol = simulate(m, [0, 0.5], 20; show_degenerated = false);
       get_alive_species(sol)
(species = [:s2], idxs = [2])

julia> sol = simulate(m, [0.5, 0], 20; show_degenerated = false);
       get_alive_species(sol)
(species = [:s1], idxs = [1])

julia> sol = simulate(m, [0.5, 0.5], 20; show_degenerated = false);
       get_alive_species(sol)
(species = [:s1, :s2], idxs = [1, 2])

julia> sol = simulate(m, [0, 0], 20; show_degenerated = false);
       get_alive_species(sol)
(species = Symbol[], idxs = Int64[])
```
"""
function get_alive_species(solution; idxs = nothing, threshold = 0)
    idxs = process_idxs(solution; idxs)
    sp = get_model(solution).species.names[idxs]
    alive = get_alive_species(solution[idxs, end]; threshold)
    (species = sp[alive], idxs = idxs[alive])
end
get_extinct_species(m::AbstractVector; threshold = 0) = findall(<=(threshold), m)
get_alive_species(m::AbstractVector; threshold = 0) = findall(>(threshold), m)
export get_alive_species, get_extinct_species
