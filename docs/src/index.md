# FiniteStateProjection.jl Documentation

```@meta
CurrentModule=FiniteStateProjection
```

## Introduction

FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). FiniteStateProjection.jl converts descriptions of reaction networks into `ODEProblem`s that can be used to compute approximate solutions of the Chemical Master Equation with packages such as [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

FiniteStateProjection.jl works by converting a `ReactionSystem` into a function that computes the right-hand side of the Chemical Master Equation:

``\frac{\mathrm{d}}{\mathrm{d} t} P(t) = A P(t)``

This function is generated dynamically via [`build_rhs`](@ref) and specialised for each `ReactionSystem`. Users can use their preferred array types supporting `CartesianIndices` and provide additional features by overloading these functions.

### Features:
- Flexible API for user-defined array types
- Automatic dimensionality reduction for systems with conserved quantities
- On-the-fly generation of specialised functions for improved performance

### Examples

#### Birth-Death Model

This example models a linear birth-death process. The reaction network is easily defined using [Catalyst.jl](https://github.com/SciML/Catalyst.jl). Our truncated state space has length 50, which is enough for this simple system.

This system has no conserved quantities, so we use a [`NaiveIndexHandler`](@ref) to map from a one-dimensional array with offset 1 to the state of the system. See [Index Handlers](@ref) for more details.

```julia
using FiniteStateProjection, DifferentialEquations

@parameters r1, r2
rs = @reaction_network begin
    r1, 0 --> A
    r2, A --> 0
end r1 r2

sys = FSPSystem(rs)

# Parameters for our system
ps = [ 10.0, 1.0 ]

# Initial values
u0 = zeros(50)
u0[1] = 1.0

prob = convert(ODEProblem, NaiveIndexHandler(sys, 1), sys, u0, 10.0, ps)
sol = solve(prob, Vern7(), atol=1e-6)
```
![Visualisation](assets/birth_death.png)

#### Telegraph Model

Here we showcase the telegraph model, a simplistic description of mRNA transcription in biological cells. We have one gene that transitions stochastically between an *on* and an *off* state and produces mRNA molecules while it is in the *on* state.

This system technically consists of three different species, namely the two states of the gene and mRNA. It is clear, however, that these are not independent as ``D_{on}(t) + D_{off}(t) = 1``. In order to solve the Chemical Master Equation we can therefore recover ``D_{off}(t)`` from the other variables and the entire state of the system is described by only two variables: ``D_{on}(t)`` and ``M(t)``, as well as the total number of genes, which is a constant equal to $1$. The default index handler class [`DefaultIndexHandler`](@ref) does this for us automatically and maps the state of the system to a two-dimensional array. This showcases that we can often reduce the number of species in the system to make it easier to solve numerically.

```julia
using FiniteStateProjection, DifferentialEquations

@parameters r1 r2 r3 r4
rs = @reaction_network begin
    r1, G_on --> G_on + M
    (r2, r3), G_on <--> G_off
    r4, M --> 0
end r1 r2 r3 r4

sys = FSPSystem(rs)

# There is one conserved quantity: G_on + G_off
cons = conservedquantities([1,0,0], sys)

# Parameters for our system
ps = [ 15.0, 0.25, 0.15, 1.0 ]

# Since G_on + G_off = const. we do not have to model the latter separately
u0 = zeros(2, 50)
u0[1,1] = 1.0

prob = convert(ODEProblem, DefaultIndexHandler(sys, 1), sys, u0, 10.0, (ps, cons))
sol = solve(prob, Vern7(), atol=1e-6)
```
![Visualisation](assets/telegraph.png)


## FSP Basics

```@docs
FSPSystem
conservationlaws
conservedquantities
```

## Index Handlers

The task of an index handler is to provide a mapping between the way the solution of the FSP is stored, usually a multidimensional array, and the states it represents. The standard approach is to store the states of a system with ``s`` reactions as an ``s``-dimensional array and have the index ``(i_1, \ldots, i_s)`` correspond to the state ``(n_1 = i_1, \ldots, n_s = i_s)``. This is implemented by the class [`NaiveIndexHandler`](@ref), which accepts an offset argument to deal with Julia's 1-based indexing (so the Julia idex $(1,\ldots,1)$ corresponds to the state with no molecules). For systems with conservation laws the [`DefaultIndexHandler`](@ref) class generally stores the data more efficiently and is the preferred choice.

User-defined index handlers should inherit from [`AbstractIndexHandler`] and implement the following methods:
- [`getsubstitutions`](@ref)
- [`build_rhs_header`](@ref)
- [`singleindices`](@ref) (optional)
- [`pairedindices`](@ref) (optional)

```@docs
AbstractIndexHandler
singleindices
pairedindices
getsubstitutions
NaiveIndexHandler
DefaultIndexHandler
reducedspecies
elidedspecies
```

## Function Building

```@docs
Base.convert
build_rhs
unpackparams
build_rhs_header
build_ratefuncs
build_rhs_firstpass
build_rhs_secondpass
```
