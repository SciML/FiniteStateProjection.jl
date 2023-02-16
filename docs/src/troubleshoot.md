# [Troubleshooting](@id troubleshoot)

Solving the Chemical Master Equation numerically is a difficult task and errors are liable to crop up. The following section presents some ways to catch common errors.

## Ensure your state space has the right dimension

If your are solving an SIR model with three species, ``S``, ``I`` and ``R``, your state space will be 3-dimensional. FiniteStateProjection.jl computes probabilities for all states simultaneously and stores the results in a 3-dimensional array. In particular, `u0` must have type `Float64` or similar as it represents numbers between 0 and 1.

```julia
# correct
u0 = zeros(101, 101, 101)
u0[100,2,1] = 1.0           # start with 1 infected and 99 susceptible individuals

# incorrect
u0 = [99,1,0]               # wrong type (Int) and shape (1D)
```

## Ensure your state space is big enough

A common reason for the solver returning nonsensical solutions is a state space that is too small. Since the Finite State Projection works with a finite-dimensional approximation of the system, the number of states considered can have a large impact on accuracy. The loss of accuracy due to using a smaller state space is the truncation error.

A good way to check whether your state space is large enough is to solve the CME until the required time ``t`` and to sum up the probabilities for each state - this gives the probability that the system will have remained in the truncated state space from start to end. If this quantity is noticeably less than 1, the state space is likely too small. As an informal rule of thumb, a value of less than 95% indicates that the solution will not be reliable.

```julia
# always true
sum(u0) == 1                # our initial condition lies in the truncated state space

# good
sum(ut) >= 0.99             # our truncation covers most of the relevant states

# bad
sum(ut) < 0.95              # we do not have enough coverage
```

The above does not work directly when computing steady-state probabilities as the value will usually drop to 0 for large enough ``t``. In this case it is a good idea to redo the computation with a larger state space - the results should agree if the truncation error is small.

## Ensure your propensities are positive

This point might seem obvious, but errors in the rate functions, or an incorrectly chosen truncation, can lead to negative reaction propensities that will typically result in numerical instabilities. As an example, consider the following version of the SI model where the population size ``S + I = N`` is constant, allowing us to rewrite the system using only one species ``I`` (with ``S = N - I``):

```julia
rn = @reaction_network begin
   σ * (N - I), I --> 2I
   ρ, I --> 0
end σ ρ N

sys_fsp = FSPSystem(rn)
```

Here the propensity function for the first reaction will negative if ``I > N``, so the following may result in numerical instabilities:

```julia
u0 = zeros(30)
u0[2] = 1

# N is too small for the state space!
prob_fsp = convert(ODEProblem, sys_fsp, u0, (0, 100.), [ 1., 1., 20 ])
```

## Ensure you are using the right solver

The Chemical Master Equation is generally very stiff and requires a solver that can handle this stiffness, see [Tips & Tricks](@ref tips). If your solver fails, first check if any of the above points apply. You may be able to get a different solver to work; this requires some experimentation. Anecdotally some systems, particularly oscillatory ones such as the Schlögl model, can pose significant challenges to most solvers and take inordinate amounts of time to solve. I am not aware of any solutions for this at the moment, but please consider opening an issue on GitHub if you encounter examples of this sort.
