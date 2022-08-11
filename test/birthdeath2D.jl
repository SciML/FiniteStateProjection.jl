using Test
using OrdinaryDiffEq
using SteadyStateDiffEq
using Distributions
using FiniteStateProjection
using SparseArrays
using LinearAlgebra
using Sundials

marg(vec; dims) = dropdims(sum(vec; dims); dims)

@parameters r1, r2, s1, s2
rs = @reaction_network begin
    r1, 0 --> A
    r2, A --> 0
    s1, 0 --> B
    s2, B --> 0
end r1 r2 s1 s2

sys = FSPSystem(rs)

prs = exp.(2 .* rand(2))
ps = [ prs[1], prs[1] / exp(4 * rand()), 
       prs[2], prs[2] / exp(4 * rand()) ]

Nmax = 130
u0 = zeros(Nmax+1, Nmax+1)
u0[1] = 1.0

tt = [ 0.25, 1.0, 10.0 ]

prob = convert(ODEProblem, sys, u0, 10.0, ps)
sol = solve(prob, Vern7(), abstol=1e-9, reltol=1e-6, saveat=tt)

@test marg(sol.u[1], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[1]))), 0:Nmax) atol=1e-4
@test marg(sol.u[1], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[1]))), 0:Nmax) atol=1e-4

@test marg(sol.u[2], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[2]))), 0:Nmax) atol=1e-4
@test marg(sol.u[2], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[2]))), 0:Nmax) atol=1e-4

@test marg(sol.u[3], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[3]))), 0:Nmax) atol=1e-4
@test marg(sol.u[3], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[3]))), 0:Nmax) atol=1e-4

A = convert(SparseMatrixCSC, sys, (Nmax+1, Nmax+1), ps, 0)
f = (du,u,t) -> mul!(du, A, u)

probA = ODEProblem(f, u0, 10.0)
solA = solve(prob, Vern7(), abstol=1e-9, reltol=1e-6, saveat=tt)

@test sol.u[1] ≈ solA.u[1] atol=1e-4
@test sol.u[2] ≈ solA.u[2] atol=1e-4
@test sol.u[3] ≈ solA.u[3] atol=1e-4

## Steady-State Tests

prob = convert(SteadyStateProblem, sys, u0, ps)
sol = solve(prob, SSRootfind())
sol.u ./= sum(sol.u)

@test marg(sol.u, dims=2) ≈ pdf.(Poisson(ps[1] / ps[2]), 0:Nmax) atol=1e-4
@test marg(sol.u, dims=1) ≈ pdf.(Poisson(ps[3] / ps[4]), 0:Nmax) atol=1e-4

A = convert(SparseMatrixCSC, sys, (Nmax+1, Nmax+1), ps, SteadyState())
f = (du,u,p,t) -> mul!(vec(du), A, vec(u))

probA = SteadyStateProblem(f, u0)
solA = solve(probA, SSRootfind())
solA.u ./= sum(solA.u)

@test sol.u ≈ solA.u atol=1e-4
