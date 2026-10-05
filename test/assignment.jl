using GraphsOptim
using LinearAlgebra
using Test

A = [7 4 5; 6 9 8; 9 9 11]
P = min_cost_assignment(-A; integer=true)
@test P ≈ I

@test GraphsOptim.is_doubly_stochastic(min_cost_assignment(rand(10, 10); integer=false))
@test GraphsOptim.is_permutation_matrix(min_cost_assignment(rand(10, 10); integer=true))

@testset "Linear sum assignment" begin
    for n in (1, 5, 20),
        C in (rand(n, n), Float64.(rand(0:2, n, n))),
        maximize in (false, true)

        σ = GraphsOptim.linear_sum_assignment(C; maximize)
        @test isperm(σ)
        P = min_cost_assignment(maximize ? -C : C; integer=true)
        @test sum(C[i, σ[i]] for i in 1:n) ≈ sum(C .* P)
    end
end
