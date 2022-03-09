using GLPK
using Graphs
using GraphsOptim
using LinearAlgebra
using SparseArrays
using Test

@testset verbose = true "GraphsOptim.jl" begin
    include("minimum_cost_flow.jl")
    include("maximum_weight_matching.jl")
end
