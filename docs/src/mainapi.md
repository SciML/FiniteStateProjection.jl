# Main API
```@meta
CurrentModule = FiniteStateProjection
```

## Reaction Systems

The following describes the [`FSPSystem`](@ref) struct, a thin wrapper around ModelingToolkit's `ReactionSystem` for use with this package. The main functionality of this type is determining conserved quantities for a reaction system.

```@docs
FSPSystem
conservationlaws(::FSPSystem)
conservedquantities
```

## Creating ODE systems

The following methods are the main way to create a system of ODEs representing the time-dependent FSP. This package provides a flexible way to represent the FSP in memory via index handlers, see [Index Handlers] for more information. 
 
```@docs
Base.convert(::Type{ODEFunction}, ::AbstractIndexHandler, ::FSPSystem)
Base.convert(::Type{ODEProblem}, ::AbstractIndexHandler, ::FSPSystem, u0, tmax, p)
```
