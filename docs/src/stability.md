# Stability

This section covers how to analyse the stability of a community once simulated
with [EcologicalNetworksDynamics](https://econetoolbox.github.io/EcologicalNetworksDynamics.jl/).
Stability can be notoriously measured in many different ways.
We provide in what follows an overview of the different stability metrics
available from our package.
We have aimed to implement the most common and used stability metrics.

To begin we need to import the companion package.

```@example econetd
using EcoNetDynOutputs # Companion package - stability functions.
```

## Simulating a simple community

First, let's consider a simple community of a consumer feeding on a producer.

```@example econetd
using EcologicalNetworksDynamics

fw = Foodweb([:consumer => :producer])
m = default_model(fw)
```

Before going further, we can check species indices, which is going to
useful for subsequent steps.

```@example econetd
m.species.index
```

We see that species 1 is the consumer, and species 2 is the producer.
We can simulate the dynamic of this simple model to find its steady-state.

```@example econetd
B0 = [1, 1] # Initial biomass.
sol = simulate(m, B0, 1_000)
Beq = sol[end]
```

At steady-state the consumer has a biomass near 0.39,
and the producer has a biomass near 0.18.

## Computing the community Jacobian

Many stability metrics are derived from the Jacobian of the dynamical system.
The Jacobian is simply a matrix whose elements specifies how species behave in face
of a small disturbance.
It tells for example if species are prone to go back quickly to their undisturbed state
or, on the contrary, if species tend to go away.

We can compute the Jacobian of the community previously simulated.
To do so, we need to specify two arguments:
- `m` specifying the model
- `B` the vector of biomass where to evaluate the jacobian
Most of the time `B` is going to be species biomass at equilibrium,
as the Jacobian is particularly useful and relevant to perform
the stability analysis near an equilibrium.

In our setting, we have

```@example econetd
j = jacobian(m, Beq)
```


Various stability metrics can be derived from the Jacobian.
The most common is by far the asymptotic resilience, formally defined as

```math
\text{resilience} = \max_i \Re(\lambda_i(J))
```


where ``\lambda_i`` denote the eigenvalues of the Jacobian ``J``.


We can compute the [`resilience`](@ref) of the system with

```@example econetd
resilience(j)
```

We expect a negative value, as we have seen that the equilibrium is stable.
The more negative the resilience, the more stable the community is.

Another common stability metric derived from the Jacobian, that is
increasingly used is the community reactivity.
While the resilience gives the long-term recovery rate of the community,
the [`reactivity`](@ref) informs on the contrary on the short-term recovery rate.
It is formally defined as

```math
\text{reactivity} = \max_i \lambda_i (\frac{J + J^T}{2})
```

It more intuitively corresponds to the degree to which a disturbance can be amplified in the worst case scenario. It can be computed with

```@example econetd
reactivity(j)
```

A positive value means that the system can go away from its equilibrium after a
disturbance, before eventually recovering.

## Sensitivity matrix

While the Jacobian describes the community recovery near its equilibrium,
after an instantaneous disturbance, it does not capture
how species respond to sustained change in environmental conditions,
such as an increase in mortality.
The response to this type of disturbances (hereafter 'press') is captured by the
*sensitivity matrix*.
Elements of the sensitivity matrix for instance quantifies
how an increase is species mortality
affects the biomass of another species.
The sensitivity matrix is simply the inverse of the interaction matrix.
Because in our model interactions are density-dependent,
the vector of species biomass (where to evaluate interactions)
should be specified.

```@example econetd
A = get_interaction(m, Beq)
```

Formally, the interaction coefficients are defined as

```math
A_{ij} = \frac{\partial (\frac{1}{B_i}\frac{\text{d} B_i}{\text{d}t})}{\partial B_j}
```

that is the partial derivative of the growth rate of species ``i``
relative to a change in biomass of species ``j``.

In our toy example, we see that species 2 (producer)
is self-regulated (`A[2, 2] < 0`).
Moreover because species 1 feed on species 2
it receives a positive interaction (`A[1, 2]`),
while species 2 receives a negative interaction (`A[2, 1]`).
The sensitivity matrix can be directly computed

The sensivitiy matrix, being the inverse of the interaction matrix,
informs on the inverse relationship

```math
S_{ij} = \frac{\partial B_i^\text{eq}}{\partial (\frac{1}{B_j}\frac{\text{d}B_j}{\text{d}t})}
```

that is how species equilibrium changes in response to variation in species growth rates. For example, how species biomass respond to an increase in mortality
in other species.


```@example econetd
S = sensitivity(m, Beq)
```

We see that an increase in the producer mortality
should not affect its biomass (`S[2, 2]=0`),
while it should decrease the biomass of the consumer (`S[1, 2]>0`).

```@example econetd
m.d[2] = 0.1 # Was set to 0 before.
Bnew = simulate(m, Beq, 1_000).u[end]
Bnew .- Beq
```

The sensitivity matrix can notably used to compute species [`resistance`](@ref).
Resistance quantifies how species respond to disturbance in others.
Formally it is defined as

```math
\text{resistance} = \frac{\Delta B_i}{\Delta d_i}
```

that is the ratio of the change in species biomass in response normalized
by its increase in mortality rate.

```@example econetd
resistance(m, Beq)
```

Here resistance is derived from the sensitivity matrix, which assumes a small disturbance so that nonlinearities can be neglected.
But resistance can be also computed from simulation with

```@example econetd
resistance_simulation(m, Beq)
```

We see first that output values are close to the one returned [`resistance`](@ref).
The slight difference is explained by the nonlinearities of the species response
which is not accounted for when using [`resistance`](@ref).
Furthermore, we see that species (that is, the consumer) is less resistant
than the plant it feeds on (species 2).
Ecologically speaking, this can be understood by the fact that
when the mortality of both species increase
the plant benefit from a release of the feeding pressure
which offset its biomass loss.
On the contrary, the consumer
suffers even more as in addition of the increase in its mortality rate
it also perceives a reduction in the amount of the plant it feeds on
as the mortality of the plant also increases.
In short, the difference in resistance value between two species translate their asymmetric relationship.


Both functions, [`resistance`](@ref), [`resistance_simulation`](@ref)
have multiple keyword arguments to specify which species is affected by the disturbance
the response of which species should be measured, and if the measure
should be aggregated (averaged) at the community level.

## Secondary extinctions and robustness

Another set of stability metrics focus on secondary extinctions (or cascading extinctions).
Secondary extinctions measure the number of extinctions
following the extinction of a first species.
It informs on the importance, or more specifically the 'keystoness', of that species.
If the extinction of that species result in many secondary extinctions,
it means that the species is key for the community (or is a keystone species)
as its presence allows many other species to coexist.
On the contrary, if the number of secondary extinctions is null it means
that the species is less important (regarding the community diversity at least).

In our toy example we can compute [`secondary_extinctions`](@ref)
as follow

```@example econetd
secondary_extinctions(m, 1, Beq) # Species 1 is set extinct.
```

```@example econetd
secondary_extinctions(m, 2, Beq) # Species 2 is set extinct.
```

We see that the extinction of the consumer leads to no secondary extinction,
while the extinction of the plant leads to the extinction of the consumer
that depends on it.

We can compute at the level of the community, the extent to which a community
is prone to secondary extinctions if a species goes extinct.
This stability measure is given by the [`robustness`](@ref).
It is defined as the inverse of the number of secondary extinction we can expect.
In our example, we have seen that we have two possibilities : 

1. The plant goes extinct first and we have one secondary extinction
2. The consumer goes extinct first and we have no secondary extinction

Therefore we expect a robustness value close to 

```math
\text{robustness}_\text{exp} = (\frac{1 + 0}{2})^{-1} = 2
```

We can check this numerically with

```@example econetd
robustness(m)
```

We do observe a value close to 2. The value is not exact because the
computation is performed on a set of random extinction sequence,
increasing this number with `n_rep` should result in the convergence
of the robustness value towards 2 in this case.
