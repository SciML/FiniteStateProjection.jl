# FiniteStateProjection.jl Documentation

```@meta
CurrentModule=FiniteStateProjection
```

## Introduction

FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). FiniteStateProjection.jl converts descriptions of reaction networks into `ODEProblem`s that can be used to compute approximate solutions of the Chemical Master Equation with packages such as [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

FiniteStateProjection.jl works by converting a `ReactionSystem` into a function that computes the right-hand side of the Chemical Master Equation:

``\frac{\mathrm{d}}{\mathrm{d} t} P(t) = A P(t)``

This function is generated dynamically as an `ODEFunction` for use with DifferentialEquations.jl and specialised for each `ReactionSystem`. Users can use their preferred array types and provide additional features by overloading the functions in this package. Alternatively the matrix `A` can be constructed as a `SparseMatrixCSC`.

## Features
- Built on top of [Catalyst.jl](https://github.com/SciML/Catalyst.jl)
- FSP equations are generated as `ODEFunction`/`ODEProblem`s and can be solved with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), with on-the-fly generation of targeted functions for improved performance
- The Chemical Master Equation can be generated as a `SparseMatrixCSC`
- Flexible API for user-defined array types via `IndexHandler`s
- Automatic dimensionality reduction for systems with conserved quantities
