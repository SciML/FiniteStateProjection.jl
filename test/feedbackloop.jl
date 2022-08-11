# Thanks to Xiaoming Fu for sharing his example code

using Test
using OrdinaryDiffEq
using SteadyStateDiffEq
using FiniteStateProjection
using SparseArrays

@parameters σ_off σ_on ρ_on ρ_off d

rs = @reaction_network begin
	σ_off, G + P → 0
	σ_on * (1 - G), 0 ⇒ G + P
	ρ_on, G → G + P
	ρ_off * (1-G), 0 ⇒ P
	d, P → 0
end σ_off σ_on ρ_on ρ_off d

ps = [ 1, 0.1, 20, 1, 1]
Nmax = 200

sys = FSPSystem(rs)

function create_A(ps, Nmax)
	σ_off, σ_on, ρ_on, ρ_off, d = ps

	A00 = sparse([ j == i ? -(ρ_off + σ_on + (i-1) * d) : j == i - 1 ? ρ_off : j == i + 1 ? i * d : 0 for i = 1:Nmax, j=1:Nmax])
	A01 = sparse([ j == i + 1 ? i * σ_off : 0 for i = 1:Nmax, j = 1:Nmax])

	A10 = sparse([ j == i - 1 ? σ_on : 0 for i = 1:Nmax, j = 1:Nmax])
	A11 = sparse([ j == i ? -(ρ_on + (i - 1) * (d + σ_off)) : j == i - 1 ? ρ_on : j == i + 1 ? i * d : 0 for i = 1:Nmax, j=1:Nmax ])

	A = SparseMatrixCSC{Float64, Int64}([ A00 A01; A10 A11 ])

	# Convert to column-major ordering used by Julia
	J = vec(reshape(1:2*Nmax, 2, Nmax)')
	P = sparse(1:2*Nmax, J, 1)
	P' * A * P
end

A = create_A(ps, Nmax)
A_fsp = convert(SparseMatrixCSC, sys, (2, Nmax), ps, 0.)

@test A ≈ A_fsp

tt = [ 1.0, 10.0, 20.0 ]

u0 = zeros(2, Nmax)
u0[1] = 1.0

prob = convert(ODEProblem, sys, u0, maximum(tt), ps)
sol = solve(prob, Vern7(), abstol=1e-9, reltol=1e-6, saveat=tt)

f = (du,u,t) -> mul!(du, A, u)

probA = ODEProblem(f, u0, 10.0)
solA = solve(prob, Vern7(), abstol=1e-9, reltol=1e-6, saveat=tt)

@test sol.u[1] ≈ solA.u[1] atol=1e-4
@test sol.u[2] ≈ solA.u[2] atol=1e-4
@test sol.u[3] ≈ solA.u[3] atol=1e-4

## Steady-State Tests

function create_A_ss(ps, Nmax)
    A = create_A(ps, Nmax)

    for i in 1:size(A, 1)
        A[i,i] -= sum(A[:,i])
    end

    A
end

A_ss = create_A_ss(ps, Nmax)
A_fsp_ss = convert(SparseMatrixCSC, sys, (2, Nmax), ps, SteadyState())

@test A_ss ≈ A_fsp_ss

prob = convert(SteadyStateProblem, sys, u0, ps)
sol = solve(prob, SSRootfind())
sol.u ./= sum(sol.u)

f = (du,u,p,t) -> mul!(vec(du), A, vec(u))

probA = SteadyStateProblem(f, u0)
solA = solve(prob, SSRootfind())
solA.u ./= sum(solA.u)

@test sol.u ≈ solA.u atol=1e-4
