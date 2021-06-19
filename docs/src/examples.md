# Examples

## Birth-Death Model

This example models a linear birth-death process. The reaction network is easily defined using [Catalyst.jl](https://github.com/SciML/Catalyst.jl). Our truncated state space has length 50, which is enough for this simple system.

```julia
using FiniteStateProjection
using OrdinaryDiffEq

@parameters σ d
rn = @reaction_network begin
    σ, 0 --> A
    d, A --> 0
end σ d

sys = FSPSystem(rn)

# Parameters for our system
ps = [ 10.0, 1.0 ]

# Initial values
u0 = zeros(50)
u0[1] = 1.0

prob = ODEProblem(sys, u0, (0, 10.0), ps)
sol = solve(prob, Vern7(), atol=1e-6)
```
![Visualisation](assets/birth_death.png)

## Telegraph Model

Here we showcase the telegraph model, a simplistic description of mRNA transcription in biological cells. We have one gene that transitions stochastically between an *on* and an *off* state and produces mRNA molecules while it is in the *on* state.

This system technically consists of three different species, namely the two states of the gene and mRNA. It is clear, however, that these are not independent as ``D_{on}(t) + D_{off}(t) = 1``. In order to solve the Chemical Master Equation we can therefore recover ``D_{off}(t)`` from the other variables and the entire state of the system is described by only two variables: ``D_{on}(t)`` and ``M(t)``, as well as the total number of genes, which is a constant equal to $1$. The index handler class [`ReducingIndexHandler`](@ref) does this for us automatically and maps the state of the system to a two-dimensional array. This showcases that we can often reduce the number of species in the system to make it easier to solve numerically.

```julia
using FiniteStateProjection
using OrdinaryDiffEq

@parameters ρ σ_on σ_off d
rn = @reaction_network begin
    ρ, G_on --> G_on + M
    (σ_on, σ_off), G_off <--> G_on
    d, M --> 0
end ρ σ_on σ_off d

# This automatically reduces the dimensionality of the
# network by exploiting conservation laws
ih = ReducedIndexHandler(rn)
sys = FSPSystem(rn, ih)

# There is one conserved quantity: G_on + G_off
cons = conservedquantities([1, 0, 0], sys)

# Parameters for our system
ps = [ 15.0, 0.25, 0.15, 1.0 ]

# In the reduced model, G_off = 1 - G_on does not have to be tracked
u0 = zeros(2, 50)
u0[1,1] = 1.0

prob = ODEProblem(sys, u0, (0, 10.0), (ps, cons))
sol = solve(prob, Vern7(), atol=1e-6)
```
![Visualisation](assets/telegraph.png)

!!! This feature is currently being moved to [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and is subject to removal from this package.
