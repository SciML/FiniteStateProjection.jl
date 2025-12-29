using Test
using FiniteStateProjection
using SparseArrays
using OrdinaryDiffEq

# Test that the core RHS evaluation is allocation-free
@testset "Allocation Tests - RHS Evaluation" begin
    # Simple 2-species model
    rs = @reaction_network begin
        k1, 0 --> A
        k2, A --> 0
        k3, A --> A + B
        k4, B --> 0
    end

    sys = FSPSystem(rs)
    pmap = [:k1 => 1.0, :k2 => 0.5, :k3 => 1.0, :k4 => 0.5]
    ps = (1.0, 0.5, 1.0, 0.5)

    Amax = 20
    Bmax = 30
    u0 = zeros(Amax + 1, Bmax + 1)
    u0[1, 1] = 1.0
    du = similar(u0)

    # Build the ODE problem to get the RHS function
    prob = ODEProblem(sys, u0, 1.0, pmap)
    f = prob.f.f

    # Warm up
    f(du, u0, ps, 0.0)

    # Test that RHS evaluation is allocation-free
    allocs = @allocated f(du, u0, ps, 0.0)
    @test allocs == 0
end

@testset "Allocation Tests - Index Iteration" begin
    ih = FiniteStateProjection.DefaultIndexHandler{3}()
    dims = (10, 20, 30)
    shift = CartesianIndex(1, 0, -1)
    u = zeros(dims)

    # Test that iteration itself is allocation-free by using a function barrier
    # This avoids measuring JIT compilation overhead
    function count_single(ih, u)
        count = 0
        for idx in FiniteStateProjection.singleindices(ih, u)
            count += 1
        end
        count
    end

    function count_paired(ih, dims, shift)
        count = 0
        for (idx_in, idx_out) in FiniteStateProjection.pairedindices(ih, dims, shift)
            count += 1
        end
        count
    end

    # Warmup
    count_single(ih, u)
    count_paired(ih, dims, shift)

    # Test singleindices is allocation-free
    allocs_single = @allocated count_single(ih, u)
    @test allocs_single == 0

    # Test pairedindices is allocation-free
    allocs_paired = @allocated count_paired(ih, dims, shift)
    @test allocs_paired == 0
end

@testset "Allocation Tests - 3-species Telegraph Model" begin
    # Full telegraph model from test suite
    rs = @reaction_network begin
        σ_on * (1 - G), 0 --> G
        σ_off, G --> 0
        ρ_m, G --> G + M
        δ_m, M --> 0
        ρ_p, M --> M + P
        δ_p, P --> 0
    end

    sys = FSPSystem(rs)
    pmap = [:σ_on => 0.5, :σ_off => 0.3, :ρ_m => 1.0, :δ_m => 0.2, :ρ_p => 2.0, :δ_p => 0.1]
    ps = (0.5, 0.3, 1.0, 0.2, 2.0, 0.1)

    Mmax = 10
    Pmax = 20
    u0 = zeros(2, Mmax + 1, Pmax + 1)
    u0[1] = 1.0
    du = similar(u0)

    prob = ODEProblem(sys, u0, 1.0, pmap)
    f = prob.f.f

    # Warm up
    f(du, u0, ps, 0.0)

    # Test allocation-free RHS
    allocs = @allocated f(du, u0, ps, 0.0)
    @test allocs == 0
end
