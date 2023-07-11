using Graphs
using GraphsOptim
using LinearAlgebra
using Test

@test GraphsOptim.is_doubly_stochastic(first(faq(rand(5, 5), rand(5, 5))))

@testset "Empty graphs" begin
    A = adjacency_matrix(SimpleGraph(5))
    B = adjacency_matrix(SimpleGraph(5))
    _, dist, converged = faq(A, B)
    @test dist ≈ 0.0
    @test converged
end

@testset "Path and complete graph" begin
    A = adjacency_matrix(path_graph(3))
    B = adjacency_matrix(complete_graph(3))
    _, dist, converged = faq(A, B)
    @test dist ≈ sqrt(2)
    @test converged
end
