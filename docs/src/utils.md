# Utilities

This module provides utility functions to simplify the post-processing and analysis of simulations from the main `EcologicalNetworksDynamics` workflow. It includes functions for extracting time series data, identifying the persistence of species, and other common tasks.

!!! note "Package Status"
    This package is under active development. Feedback and requests for new utility functions are highly welcome. Please open an issue on our GitHub repository to suggest improvements.

```julia
# Loading packages
using EcologicalNetworksDynamics, EcoNetPostProcessing
using StatsBase
```

## Extract Time Series

A common task is to analyze the population dynamics of species over time, particularly focusing on the steady state reached at the end of a simulation. The [`extract_last_timesteps`](@ref) function is the primary tool for this.

### `extract_last_timesteps`

```@docs
extract_last_timesteps
```

**Example: Analyzing Final Community State**

This example shows how to extract the last 100 time steps of a simulation to analyze the stable state of the community.

```julia
# Simulate a simple food web
fw = Foodweb([0 0; 1 0]) # A -> B
m = default_model(fw)
B0 = [1.0, 0.5] # Initial biomass
t_end = 10_000
sol = simulate(m, B0, t_end)

# Extract the last 100 timesteps for all species
final_state = extract_last_timesteps(sol; last = 100)

# Calculate the mean biomass of each species at equilibrium
mean_equilibrium_biomass = mean(final_state, dims=2)

# Extract only the last 10% of the simulation for species "s2" (the consumer)
consumer_final = extract_last_timesteps(sol; last = "10%", idxs = "s2")
```

**Key Points:**

* Use the `last` argument to specify either a fixed number of time steps (`100`) or a percentage of the total simulation (`"10%"`).
* The `idxs` argument allows you to select specific species by their index (`[2, 1]`) or name (`["s2", "s1"]`).
* The function returns a matrix where rows are species and columns are time points.

---

## Identify Persisting and Extinct Species

After a simulation, it's crucial to determine which species survived and which went extinct. The following functions help you identify species based on a biomass threshold at the end of the simulation.

### `get_alive_species`

```@docs
get_alive_species
```

### `get_extinct_species`

While not directly exported for a solution (see note below), the logic for a single biomass vector is available. For a full solution, you can use the output of `get_alive_species` to infer extinct species.

**Example: Classifying Species by their Fate**

```julia
fw = Foodweb([0 0 0; 1 0 0; 0 1 0]) # A -> B -> C
m = default_model(fw)
sol = simulate(m, [1.0, 0.0, 0.2], 1_000; extinction_threshold = 1e-5) # B begins at 0, so B its consumer C will go extinct
# Get the alive species at the end
alive_species = get_alive_species(sol)
println("Alive species: $(alive_species.species) at indices $(alive_species.idxs)")
# To find extinct species, we can get all species and find the difference
all_species_indices = 1:size(sol, 1) # Indices 1 to N
extinct_species_indices = setdiff(all_species_indices, alive_species.idxs)
println("Extinct species indices: $extinct_species_indices")
# You can also check the final biomass directly
final_biomass = sol.u[end]
println("Final biomass: $final_biomass")
# Species with biomass <= 0 are extinct.
```

```julia
fw = Foodweb([0 0 0; 1 0 0; 0 1 0]) # A -> B -> C
m = default_model(fw)
sol = simulate(m, [1.0, 0.0, 0.2], 1_000) # B begins at 0, so B its consumer C will go extinct
# Get the alive species at the end, but Species 3 is not considered extinct
# Because we did not set extinction threshold in simulate() and in get_alive_species()
alive_species = get_alive_species(sol)
println("Alive species: $(alive_species.species) at indices $(alive_species.idxs)")
# C will be also considered extinct if we set an extinction threshold
alive_species = get_alive_species(sol; threshold = 1e-5)
println("Alive species: $(alive_species.species) at indices $(alive_species.idxs)")

# Also work with a vector of biomass:.
get_alive_species(vec(extract_last_timesteps(sol; last = 1)))
# You can also use the reciprocal function ``get_extinct_species` on biomass vector:
get_extinct_species(vec(extract_last_timesteps(sol; last = 1)))
```

**Key Points:**

* The `threshold` argument allows you to define the biomass level below which a species is considered extinct (default is `0`).
* `get_alive_species` returns a named tuple `(species, idxs)` for easy access to both the names and the indices of surviving species.
* For a simple vector of biomass values, you can use `get_alive_species(vector)` and `get_extinct_species(vector)` directly.

---

## Other Utility Functions

These functions are used internally by the main utilities but can also be helpful for advanced users building custom analysis pipelines.

### `process_idxs`

This function sanitizes and validates user-provided species indices or names, converting them into a consistent format for internal use.

```@docs
process_idxs
```

### `get_extinction_timesteps`

For more detailed extinction analysis, this function identifies not just *if* a species went extinct, but *when* it happened during the simulation.

**Example: Find When Extinctions Occurred**

```julia
# after a simulation
fw = Foodweb([0 0 0; 1 0 0; 0 1 0]) # A -> B -> C
m = default_model(fw)
sol = simulate(m, [1.0, 0.0, 0.2], 1_000, extinction_threshold = 1e-5)
extinction_data = EcoNetPostProcessing.get_extinction_timesteps(sol)
for (i, sp) in enumerate(extinction_data.species)
    println("Species $sp (index $(extinction_data.idxs[i])) went extinct at timestep $(extinction_data.extinction_timestep[i])")
end
```
