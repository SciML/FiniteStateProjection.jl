# FiniteStateProjection.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kaandocal.github.io/FiniteStateProjection.jl/dev/)

Finite State Projection [[1]](#1)  algorithms for chemical reaction networks based on [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl). Converts descriptions of reaction networks into `ODEProblem`s for use with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

## Features:
- Automatic dimensionality reduction for systems with conserved quantities
- On-the-fly generation of targeted functions for improved performance

## Examples
### Birth-Death System
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
![Visualisation](docs/src/assets/birth_death.png)

### Telegraph Model
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
![Visualisation](docs/src/assets/telegraph.png)

## TODO:
- Add bursty reactions
- Add stationary FSP support

## References

<a id="1">[1]</a> B. Munsky and M. Khammash, "The Finite State Projection algorithm for the solution of the Chemical Master Equation", Journal of Chemical Physics 124, 044104 (2006). https://doi.org/10.1063/1.2145882
