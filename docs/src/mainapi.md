# Main API

```@meta
CurrentModule = FiniteStateProjection
```

## Reaction Systems

This section describes the [`FSPSystem`](@ref) struct, a thin wrapper around Catalyst's `ReactionSystem` for use with this package.

```@docs
FSPSystem
```

## Creating ODE systems

The following methods convert a reaction network into a system of ODEs representing the time-dependent FSP. This package provides a flexible way to represent the FSP in memory via index handlers, see [Index Handlers] for more information.

```@docs
Base.convert(::Type{SciMLBase.ODEFunction}, ::FSPSystem)
SciMLBase.ODEProblem(::FSPSystem, u0, tmax, p)
```

## Steady-State Problems

Computing steady-state distributions can be done using the SteadyStateDiffEq.jl package. At the moment FiniteStateProjection.jl adjusts the rate matrix so that reactions leaving the truncated state space have propensity 0. Pass a [`SteadyState`](@ref) instance to the conversion methods to request this steady-state formulation.

```@docs
SteadyState
Base.convert(::Type{ODEFunction}, ::FSPSystem, ::SteadyState)
SciMLBase.SteadyStateProblem(::FSPSystem, u0, p)
```

## Reexported model-definition and problem interface

The documented workflow starts from a [Catalyst](https://docs.sciml.ai/Catalyst/stable/)
reaction network and ends in a SciML problem:

```julia
using FiniteStateProjection
using OrdinaryDiffEq

rn = @reaction_network begin
    σ, 0 --> A
    d, A --> 0
end

sys = FSPSystem(rn)
prob = ODEProblem(sys, u0, (0, 10.0), ps)
sol = solve(prob, Vern7())
```

`using FiniteStateProjection` brings the names that workflow touches into scope.
FiniteStateProjection only reexports them — they are owned and documented
upstream, at the links below.

  - Defining and inspecting the reaction network, owned by
    [Catalyst](https://docs.sciml.ai/Catalyst/stable/): `@reaction_network`,
    `ReactionSystem` (the type [`FSPSystem`](@ref) wraps), and the accessors
    `species`, `numspecies` and `reactions`, which give the number and order of
    the state-space dimensions a truncation has to cover
  - Building the problem, owned by
    [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/) and
    [CommonSolve](https://github.com/SciML/CommonSolve.jl): `ODEProblem`,
    `SteadyStateProblem`, `ODEFunction`, `solve`. FiniteStateProjection defines
    the `FSPSystem` methods of the first three, documented above.

Anything else from Catalyst must be imported from Catalyst directly. In
particular, FiniteStateProjection does not reexport Catalyst's symbolic
model-building surface (`@species`, `@parameters`, `@variables`, `@reaction`,
`Reaction`, network-composition and analysis functions), nor the other problem
types Catalyst exposes (`SDEProblem`, `JumpProblem`, `DiscreteProblem`), none of
which the finite state projection uses. Solver algorithms are not reexported
either — `solve(prob, Vern7())` still needs
[OrdinaryDiffEq](https://docs.sciml.ai/DiffEqDocs/stable/), and
`SparseMatrixCSC` for the [matrix conversions](@ref matrix_conversions) comes
from the `SparseArrays` standard library.
