using Test
using Catalyst: @reaction_network
using FiniteStateProjection: AbstractIndexHandler, DefaultIndexHandler, FSPSystem,
    PermutingIndexHandler
import FiniteStateProjection:
    build_rhs_header, getsubstitutions, pairedindices, singleindices
using SciMLBase: ODEProblem
import Base: LinearIndices, vec
import SparseArrays: SparseMatrixCSC

struct DelegatingIndexHandler{N} <: AbstractIndexHandler
    delegate::DefaultIndexHandler{N}
end

DelegatingIndexHandler{N}() where {N} = DelegatingIndexHandler(DefaultIndexHandler{N}())

singleindices(ih::DelegatingIndexHandler, arr::AbstractArray) = singleindices(ih.delegate, arr)
singleindices(ih::DelegatingIndexHandler, dims::Tuple) = singleindices(ih.delegate, dims)
pairedindices(ih::DelegatingIndexHandler, arr::AbstractArray, shift::CartesianIndex) =
    pairedindices(ih.delegate, arr, shift)
pairedindices(ih::DelegatingIndexHandler, dims::Tuple, shift::CartesianIndex) =
    pairedindices(ih.delegate, dims, shift)
getsubstitutions(ih::DelegatingIndexHandler, rs; state_sym) =
    getsubstitutions(ih.delegate, rs; state_sym)
LinearIndices(ih::DelegatingIndexHandler, arr) = LinearIndices(ih.delegate, arr)
vec(ih::DelegatingIndexHandler, arr) = vec(ih.delegate, arr)

function build_rhs_header(::FSPSystem{<:DelegatingIndexHandler})
    return quote
        birth, death = p
    end
end

@testset "Public index-handler interface" begin
    rs = @reaction_network begin
        birth, 0 --> X
        death, X --> 0
    end
    pmap = [:birth => 2.0, :death => 1.0]
    u0 = [1.0, 0.0, 0.0]

    default_system = FSPSystem(rs)
    custom_system = FSPSystem(rs, DelegatingIndexHandler{1}())

    default_matrix = SparseMatrixCSC(default_system, size(u0), pmap, 0.0)
    custom_matrix = SparseMatrixCSC(custom_system, size(u0), pmap, 0.0)
    @test custom_matrix == default_matrix

    default_problem = ODEProblem(default_system, u0, (0.0, 1.0), pmap)
    custom_problem = ODEProblem(custom_system, u0, (0.0, 1.0), pmap)
    default_du = similar(u0)
    custom_du = similar(u0)
    default_problem.f(default_du, u0, default_problem.p, 0.0)
    custom_problem.f(custom_du, u0, custom_problem.p, 0.0)
    @test custom_du == default_du
end

@testset "Index handler error handling" begin
    rs = @reaction_network begin
        k1, 0 --> A
        k2, A --> 0
        k3, 0 --> B
    end

    @testset "PermutingIndexHandler rejects a wrong number of species" begin
        @test_throws DimensionMismatch PermutingIndexHandler(rs, [:A])
    end

    @testset "PermutingIndexHandler rejects an unknown species" begin
        @test_throws ErrorException PermutingIndexHandler(rs, [:A, :C])
    end

    @testset "PermutingIndexHandler rejects a repeated species" begin
        @test_throws ArgumentError PermutingIndexHandler(rs, [:A, :A])
    end

    @testset "pairedindices rejects a dimension tuple of the wrong length" begin
        ih = DefaultIndexHandler{2}()
        @test_throws DimensionMismatch pairedindices(ih, (3,), CartesianIndex(1, 0))
        @test_throws DimensionMismatch pairedindices(ih, (3, 3, 3), CartesianIndex(1, 0))
    end

    @testset "pairedindices accepts a dimension tuple of the right length" begin
        ih = DefaultIndexHandler{2}()
        @test length(collect(pairedindices(ih, (3, 3), CartesianIndex(1, 0)))) == 6
    end
end
