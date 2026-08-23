using SciMLTesting, FiniteStateProjection, Test

# The model-definition and problem interface FiniteStateProjection deliberately
# reexports so that `using FiniteStateProjection` on its own is enough to write the
# documented workflow: build a reaction network with `@reaction_network`, wrap it in
# an `FSPSystem`, and turn that into an `ODEProblem` or `SteadyStateProblem`. Owned
# and documented upstream; kept in sync with the reexport `export` block in
# src/FiniteStateProjection.jl and the "Reexported model-definition and problem
# interface" section of docs/src/mainapi.md.
const REEXPORTS = (
    # Defining and inspecting the reaction network (Catalyst).
    Symbol("@reaction_network"), :ReactionSystem, :numspecies, :reactions, :species,
    # Building the problem; FiniteStateProjection defines the `FSPSystem` methods of
    # `ODEProblem`, `SteadyStateProblem` and the `ODEFunction` conversion.
    :ODEFunction, :ODEProblem, :SteadyStateProblem, :solve,
)

run_qa(FiniteStateProjection; reexports_allow = REEXPORTS)

@testset "Reexport surface" begin
    # Every approved reexport must actually be reachable from `using
    # FiniteStateProjection`, so the allow-list cannot drift into approving names the
    # package no longer provides. `isdefined(@__MODULE__, ...)` tests the property
    # directly: this file's `using FiniteStateProjection` is what has to bring the name
    # into scope.
    @testset "$name" for name in REEXPORTS
        @test name in names(FiniteStateProjection)
        @test isdefined(@__MODULE__, name)
    end
end
