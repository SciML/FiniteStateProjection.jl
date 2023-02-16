# [Tips and Tricks](@id tips)

The FSP approximates an infinite-dimensional system of equations by truncating it to a finite number of variables. The accuracy of the FSP therefore depends on how many variables are retained, ie.~what portion of the state space is modelled. While simple reaction networks with 1 or 2 species are not too difficult to handle using the FSP, a naive approach will require unfeasibly large truncations for systems of even moderate complexity.

## Solving Linear Equations

The Chemical Master Equation is a set of linear ODEs, and there is a vast literature on efficiently solving these. FiniteStateProjection.jl currently offers two main approaches:
* Convert a reaction system into an `ODEFunction` for use with DifferentialEquations.jl
* Convert a reaction system into a sparse matrix for use with any linear ODE solver

While there are no official benchmarks yet the first method is likely going to be more efficient for systems with time-dependent reaction rates, while the second is recommended for time-homogeneous systems.

## Picking a Solver

The FSP generally results in [stiff equations](https://en.wikipedia.org/wiki/Stiff_equation), so using a stiff solver in generally is recommended. At the time of writing it seems like the best all-round solver is `CVODE_BDF` with the `GMRES` linear solver (requires [Sundials.jl](https://github.com/SciML/Sundials.jl)):
```julia
using Sundials

(...)

sol = solve(prob, CVODE_BDF(linear_solver=:GMRES), atol=1e-6) # very fast
```

For time-homogeneous systems where the right-hand side of the CME does not depend on `t`, solving the FSP is equivalent to computing the exponential-vector product `exp(A .* t) * vec(u0)`, where `A` is the evolution operator. Here we use the `vec` function to flatten the state space into a vector that can be multiplied by the matrix `A`. We can compute the solution very efficiently using the `expmv` function provided by [Expokit.jl](https://github.com/acroy/Expokit.jl):
```julia
using Expokit

(...)

ut = similar(u0)
expmv!(vec(ut), t, A, vec(u0), tol=1e-6)  # really fast
```

## Further Comments

Choosing the right solver for large systems of ODEs can result in time savings on the order of 10-100x, and it is recommended that you experiment with a few solvers to see which works best in your case. This section is still work in progress and there has been a lot of research on accelerating the FSP and extending it to larger reaction networks which will hopefully be reviewed here soon. Feel free to share any comments or suggestions in this direction at the [Github repository](https://github.com/kaandocal/FiniteStateProjection.jl)!
