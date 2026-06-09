using FiniteStateProjection
using Aqua
using JET
using Test

@testset "Aqua" begin
    Aqua.test_all(FiniteStateProjection)
end

@testset "JET" begin
    JET.test_package(FiniteStateProjection; target_defined_modules = true)
end
