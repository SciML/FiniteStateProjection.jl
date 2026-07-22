using FiniteStateProjection
using Catalyst
using SciMLBase: ODEProblem
using SparseArrays: SparseMatrixCSC
using Symbolics: @variables
using Test
using TOML

@testset "Symbolic API and compatibility floors" begin
    @variables p

    @test FiniteStateProjection._parameter_symbol(p) == :p

    project = TOML.parsefile(joinpath(pkgdir(FiniteStateProjection), "Project.toml"))
    compat = project["compat"]
    @test compat["Catalyst"] == "16.2"
    @test compat["ModelingToolkitBase"] == "1.17"
    @test compat["SciMLBase"] == "2.144, 3"
    @test compat["SymbolicIndexingInterface"] == "0.3.43"
    @test compat["Symbolics"] == "7.13"
    @test compat["SymbolicUtils"] == "4.18"

    rs = @reaction_network begin
        birth, 0 --> A
        death, A --> 0
    end
    sys = FSPSystem(rs)
    @test all(
        ex -> ex isa Expr,
        FiniteStateProjection.build_ratefuncs(rs, sys.ih; state_sym = :idx_in)
    )
    @test sys.rfs[1]([2], 0.0, 3.0, 4.0) == 3.0
    @test sys.rfs[2]([2], 0.0, 3.0, 4.0) == 4.0

    pmap = [:birth => 3.0, :death => 4.0]
    u = [0.25, 0.5, 0.25]
    du = similar(u)
    prob = ODEProblem(sys, u, (0.0, 1.0), pmap)
    prob.f(du, u, last.(pmap), 0.0)
    A = SparseMatrixCSC(sys, (length(u),), pmap, 0)
    @test du == A * u
end
