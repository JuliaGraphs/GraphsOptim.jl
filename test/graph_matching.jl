using Graphs
using GraphsOptim
using LinearAlgebra
using Test

@testset verbose = true "FAQ" begin
    @test GraphsOptim.is_permutation_matrix(
        first(graph_matching(FAQ(), rand(5, 5), rand(5, 5)))
    )

    @testset "Empty graphs" begin
        A = adjacency_matrix(SimpleGraph(5))
        B = adjacency_matrix(SimpleGraph(5))
        _, dist, converged = graph_matching(FAQ(), A, B)
        @test dist ≈ 0.0
        @test converged
    end

    @testset "Path and complete graph" begin
        A = adjacency_matrix(path_graph(3))
        B = adjacency_matrix(complete_graph(3))
        _, dist, converged = graph_matching(FAQ(), A, B)
        @test dist ≈ sqrt(2)
        @test converged
    end
end

@testset verbose = true "GOAT" begin
    @test GraphsOptim.is_permutation_matrix(
        first(graph_matching(GOAT(), rand(5, 5), rand(5, 5)))
    )

    @testset "Path and complete graph" begin
        A = adjacency_matrix(path_graph(3))
        B = adjacency_matrix(complete_graph(3))
        _, dist, converged = graph_matching(GOAT(), A, B)
        @test dist ≈ sqrt(2)
        @test converged
    end
end
