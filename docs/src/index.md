# FiniteStateProjection.jl Documentation

```@meta
CurrentModule=FiniteStateProjection
```

## Introduction

FiniteStateProjection.jl is a package that implements Finite State Projection algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). FiniteStateProjection.jl converts descriptions of reaction networks into `ODEProblem`s that can be used to compute approximate solutions of the Chemical Master Equation with packages such as [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

The Chemical Master Equation allows us to compute the probability of a reaction system to be in a given state ``\vec n`` at a time ``t``. The state of a system with ``s`` species is described by a vector ``\vec n = (n_1, n_2, \ldots, n_s)`` of length ``s``, where ``n_1, n_2, \ldots`` are the abundances of the different species. The reaction system itself can be represented by an ``s``-dimensional array ``P[n_1,\ldots,n_s]`` that describes the probability of the system being in state ``(n_1,\ldots,n_s)``. Here ``P`` depends on time as the system evolves. Given initial conditions ``P(0)``, the time evolution of ``P`` is given by the Chemical Master Equation, which is a set of linear ODEs that can be solved numerically to predict the probability distribution of the system in the future. 

Since the number of states for a system can be infinite, the Finite State Projection approximates the Chemical Master Equation by only considering a finite number of states, namely those which are most probable. There are various ways to truncate the state space, but the most common is to provide an upper limit for each species, that is, to require ``n_1 < M_1``, ``n_2 < M_2``, etc. for fixed thresholds ``M_1, M_2, \ldots``. The number of states considered will be ``M_1 \times M_2 \times \ldots \times M_s``, and ``P`` can be represented as an array with those dimensions. While truncating the CME allows us to solve the equations numerically, the amount of space required to store ``P`` increases quickly as more species are added. 

FiniteStateProjection.jl works by converting a `ReactionSystem` into a function that computes the right-hand side of the (truncated) Chemical Master Equation:

``\frac{\mathrm{d}}{\mathrm{d} t} P(t) = A P(t)``

This function is generated dynamically as an `ODEFunction` for use with DifferentialEquations.jl and specialised for each `ReactionSystem`. Users can use their preferred array types and provide additional features by overloading the functions in this package. Alternatively the matrix `A` can be constructed as a `SparseMatrixCSC`. FiniteStateProjection.jl supports both the time-dependent Chemical Master Equation and its steady-state version.

## Features
- Built on top of [Catalyst.jl](https://github.com/SciML/Catalyst.jl)
- FSP equations are generated as `ODEFunction`/`ODEProblem`s and can be solved with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), with on-the-fly generation of targeted functions for improved performance
- The Chemical Master Equation can be generated as a `SparseMatrixCSC`
- Flexible API for user-defined array types via `IndexHandler`s

## Acknowledgments

Special thanks to [Xiamong Fu](https://github.com/palmtree2013), [Brian Munsky](https://www.engr.colostate.edu/~munsky/) and [Huy Vo](https://github.com/voduchuy) for their examples and suggestions and for contributing most of the content in the [Tips & Tricks](@ref tips) page!
