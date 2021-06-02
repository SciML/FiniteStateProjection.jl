using Test
using OrdinaryDiffEq
using Distributions
using FiniteStateProjection
using SparseArrays

marg(vec; dims) = dropdims(sum(vec; dims); dims)

@time begin
    @testset "BirthDeath2D" begin
        @parameters r1, r2, s1, s2
        rs = @reaction_network begin
            r1, 0 --> A
            r2, A --> 0
            s1, 0 --> B
            s2, B --> 0
        end r1 r2 s1 s2

        sys = FSPSystem(rs)
        
        prs = exp.(2 .* rand(2))
        ps = [ prs; prs ./ exp.(4 .* rand(2)) ]
                
        Nmax = 250
        u0 = zeros(Nmax+1, Nmax+1)
        u0[1] = 1.0
                
        tt = [ 0.25, 1.0, 10.0 ]

        prob = convert(ODEProblem, NaiveIndexHandler(sys, 1), sys, u0, 10.0, ps)
        sol = solve(prob, Vern7(), atol=1e-6, saveat=tt)

        @test marg(sol.u[1], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[1]))), 0:Nmax) atol=1e-4
        @test marg(sol.u[1], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[1]))), 0:Nmax) atol=1e-4

        @test marg(sol.u[2], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[2]))), 0:Nmax) atol=1e-4
        @test marg(sol.u[2], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[2]))), 0:Nmax) atol=1e-4

        @test marg(sol.u[3], dims=2) ≈ pdf.(Poisson(ps[1] / ps[2] * (1 - exp(-ps[2] * tt[3]))), 0:Nmax) atol=1e-4
        @test marg(sol.u[3], dims=1) ≈ pdf.(Poisson(ps[3] / ps[4] * (1 - exp(-ps[4] * tt[3]))), 0:Nmax) atol=1e-4
        
        A = convert(SparseMatrixCSC, NaiveIndexHandler(sys, 1), sys, (Nmax+1, Nmax+1), ps)
        f = (du,u,t) -> mul!(du, A, u)
        
        probA = ODEProblem(f, u0, 10.0)
        solA = solve(prob, Vern7(), atol=1e-6, saveat=tt)
        
        @test sol.u[1] ≈ solA.u[1] atol=1e-4
        @test sol.u[2] ≈ solA.u[2] atol=1e-4
        @test sol.u[3] ≈ solA.u[3] atol=1e-4
    end
    
    # Thanks to Xiaoming Fu for sharing his example code
    @testset "FeedbackLoop" begin
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
        B = similar(A)
        fill!(B, NaiveIndexHandler(1), sys, (2, Nmax), ps)
        
        A_fsp = convert(SparseMatrixCSC, NaiveIndexHandler(1), sys, (2, Nmax), ps)
        
        @test A ≈ A_fsp
        @test B ≈ A_fsp
        
        tt = [ 1.0, 10.0, 20.0 ]
        
        u0 = zeros(2, Nmax)
        u0[1] = 1.0

        prob = convert(ODEProblem, NaiveIndexHandler(sys, 1), sys, u0, maximum(tt), ps)
        sol = solve(prob, Vern7(), atol=1e-6, saveat=tt)

        f = (du,u,t) -> mul!(du, A, u)
        
        probA = ODEProblem(f, u0, 10.0)
        solA = solve(prob, Vern7(), atol=1e-6, saveat=tt)
        
        @test sol.u[1] ≈ solA.u[1] atol=1e-4
        @test sol.u[2] ≈ solA.u[2] atol=1e-4
        @test sol.u[3] ≈ solA.u[3] atol=1e-4
    end
end
