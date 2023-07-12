using Graphs
using GraphsOptim
using LinearAlgebra
using Test

@testset "FAQ" begin
    @test GraphsOptim.is_permutation_matrix(first(faq(rand(5, 5), rand(5, 5))))

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

end

@testset "GOAT" begin
    @test GraphsOptim.is_permutation_matrix(first(goat(rand(5, 5), rand(5, 5))))

    @testset "Empty graphs" begin
        A = adjacency_matrix(SimpleGraph(5))
        B = adjacency_matrix(SimpleGraph(5))
        _, dist, converged = goat(A, B)
        @test dist ≈ 0.0
        @test converged
    end


    @testset "Path and complete graph" begin
        A = adjacency_matrix(path_graph(3))
        B = adjacency_matrix(complete_graph(3))
        _, dist, converged = goat(A, B)
        @test dist ≈ sqrt(2)
        @test converged
    end

end