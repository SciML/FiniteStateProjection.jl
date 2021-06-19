# FiniteStateProjection.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaandocal.github.io/FiniteStateProjection.jl/dev/)

Finite State Projection [[1]](#1)  algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). Converts descriptions of reaction networks into `ODEProblem`s for use with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Features:
- Built on top of [Catalyst.jl](https://github.com/SciML/Catalyst.jl)
- FSP equations are generated as `ODEFunction`/`ODEProblem`s and can be solved with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), with on-the-fly generation of targeted functions for improved performance
- The Chemical Master Equation can be represented as a `SparseMatrixCSC`
- Automatic dimensionality reduction for systems with conserved quantities

More information is available in the [documentation](https://kaandocal.github.io/FiniteStateProjection.jl/dev/). Please feel free to open issues and submit pull requests! 

## Examples
### Birth-Death System
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
![Visualisation](docs/src/assets/birth_death.png)

### Telegraph Model
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
![Visualisation](docs/src/assets/telegraph.png)

## TODO:
- Add bursty reactions
- Add stationary FSP support
- Add support for sparse Jacobians

## References

<a id="1">[1]</a> B. Munsky and M. Khammash, "The Finite State Projection algorithm for the solution of the Chemical Master Equation", Journal of Chemical Physics 124, 044104 (2006). https://doi.org/10.1063/1.2145882
