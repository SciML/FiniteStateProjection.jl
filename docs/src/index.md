# FiniteStateProjection.jl Documentation

```@meta
CurrentModule=FiniteStateProjection
```

## Introduction

FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). FiniteStateProjection.jl converts descriptions of reaction networks into `ODEProblem`s that can be used to compute approximate solutions of the Chemical Master Equation with packages such as [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

FiniteStateProjection.jl works by converting a `ReactionSystem` into a function that computes the right-hand side of the Chemical Master Equation:

``\frac{\mathrm{d}}{\mathrm{d} t} P(t) = A P(t)``

This function is generated dynamically as an `ODEFunction` for use with DifferentialEquations.jl and specialised for each `ReactionSystem`. Users can use their preferred array types and provide additional features by overloading the functions in this package.

## Features
- Flexible API for user-defined array types via `IndexHandler`s
- Automatic dimensionality reduction for systems with conserved quantities
- On-the-fly generation of specialised functions for improved performance
