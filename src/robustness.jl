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

function robustness(m::Model; t_end = 1_000, n_rep = 100, threshold = 1e-6)
    S = m.richness
    B0 = simulate(m, fill(1, S), t_end; show_degenerated = false).u[end]
    rob_vec = Any[undef for _ in 1:n_rep]
    for k in 1:n_rep
        @info k
        B0_copy = deepcopy(B0)
        seq = shuffle(1:S)
        nb_ext_list = []
        i = 1
        while any(B0_copy .> 0)
            extinct_sp = seq[i]
            B0_copy[extinct_sp] = 0
            sol = simulate(m, B0_copy, t_end; show_degenerated = false)
            B0_copy = sol.u[end]
            B0_copy[B0_copy.<threshold] .= 0
            sec_ext = setdiff(findall(==(0), B0_copy), seq[1:i])
            nb_ext = length(sec_ext)
            push!(nb_ext_list, nb_ext)
            i += 1
        end
        rob_vec[k] = mean(nb_ext_list)
    end
    mean(rob_vec)
end
export robustness

