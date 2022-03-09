using GLPK
using Graphs
using GraphsOptim
using LinearAlgebra
using SparseArrays
using Test

@testset verbose = true "GraphsOptim.jl" begin
    include("mincost_flow.jl")
end
