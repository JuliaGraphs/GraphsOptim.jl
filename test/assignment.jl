using GraphsOptim
using LinearAlgebra
using Test

A = [7 4 5; 6 9 8; 9 9 11]
P = min_cost_assignment(-A; integer=true)
@test P ≈ I

@test GraphsOptim.is_doubly_stochastic(min_cost_assignment(rand(10, 10); integer=false))
@test GraphsOptim.is_permutation_matrix(min_cost_assignment(rand(10, 10); integer=true))
