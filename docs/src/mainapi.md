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
Base.convert(::Type{ODEFunction}, ::FSPSystem)
DiffEqBase.ODEProblem(::FSPSystem, u0, tmax, p)
```

## Steady-State Problems

Computing steady-state distributions can be done using the SteadyStateDiffEq.jl package. At the moment FiniteStateProjection.jl adjusts the rate matrix so that reactions leaving the truncated state space have propensity 0.

```@docs
Base.convert(::Type{ODEFunction}, ::FSPSystem, ::SteadyState)
DiffEqBase.SteadyStateProblem(::FSPSystem, u0, p)
```
