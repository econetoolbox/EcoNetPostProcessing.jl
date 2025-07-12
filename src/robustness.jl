"""
    secondary_extinctions(
    m::Model,
    extinct_sp::Integer,
    B_start::AbstractVector;
    t_end = 1_000,
    threshold = 1e-6,

)

Compute secondary extinctions following the primary extinction
of `extinct_sp` of a model `m` with initial biomasses `B_start`.

### Keyword arguments

  - `t_end` specifies the duration of simulation
  - `threshold` gives the biomass below which a species is considered extinct.

See also .
"""
function secondary_extinctions(
    m::Model,
    extinct_sp::Integer,
    B_start::AbstractVector;
    t_end = 1_000,
    threshold = 1e-6,
)
    B_copy = deepcopy(B_start)
    B_copy[extinct_sp] = 0
    sol = simulate(m, B_copy, t_end; show_degenerated = false)
    B_end = sol.u[end]
    setdiff(findall(<=(threshold), B_end), extinct_sp)
end
export secondary_extinctions

