using BasicBVP1D
using Test
using LinearAlgebra

@testset "Dense" begin
    @test include("dense_test1.jl")
    @test include("dense_test2.jl")
    @test include("dense_test3.jl")
end
@testset "Sparse" begin
    @test include("sparse_test1.jl")
    @test include("sparse_test2.jl")
    @test include("sparse_test3.jl")
end