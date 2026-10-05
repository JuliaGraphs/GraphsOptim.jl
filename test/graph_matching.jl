using Graphs
using GraphsOptim
using HiGHS
using LinearAlgebra
using SparseArrays
using Test

@testset "Step size" begin
    n = 6
    A, B = rand(n, n), rand(n, n)
    P = GraphsOptim.flat_doubly_stochastic(n)
    Q = Matrix(1.0I, n, n)[[2, 4, 6, 1, 3, 5], :]
    f(α) = tr(A * (α * P + (1 - α) * Q) * B' * (α * P + (1 - α) * Q)')
    α = GraphsOptim.graph_matching_step_size(A, B, P, Q)
    @test 0 <= α <= 1
    @test f(α) >= maximum(f, range(0, 1; length=1001)) - 1e-10
end

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

@testset "Sparse, dense and JuMP assignment" begin
    A, B = sprand(30, 30, 0.2), sprand(30, 30, 0.2)
    for algo in (FAQ(), GOAT())
        P_sparse, dist_sparse, _ = graph_matching(algo, A, B)
        P_dense, dist_dense, _ = graph_matching(algo, Matrix(A), Matrix(B))
        @test P_sparse == P_dense
        @test dist_sparse ≈ dist_dense ≈ norm(A * P_dense - P_dense * B)
    end
    P_jump, _, _ = graph_matching(FAQ(), A, B; optimizer=HiGHS.Optimizer)
    @test GraphsOptim.is_permutation_matrix(P_jump)
end
