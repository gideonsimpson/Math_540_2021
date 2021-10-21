using Poisson2D
using Test
using LinearAlgebra

@testset "Direct Solve" begin
    @test include("direct_test1.jl")
end
