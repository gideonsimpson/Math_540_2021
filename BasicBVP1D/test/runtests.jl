using BasicBVP1D
using Test
using LinearAlgebra

@testset "Dense" begin
    @test include("dense_test1.jl")
    @test include("dense_test2.jl")
end