using Aqua
using Documenter
using GLPK
using Graphs
using GraphsOptim
using JuliaFormatter
using LinearAlgebra
using SparseArrays
using Test

@testset verbose = true "GraphsOptim.jl" begin
    @testset verbose = true "Code quality (Aqua.jl)" begin
        Aqua.test_all(GraphsOptim; ambiguities=false)
    end

    @testset verbose = true "Code formatting (JuliaFormatter.jl)" begin
        @test format(GraphsOptim; verbose=false, overwrite=false)
    end

    @testset verbose = true "Doctests (Documenter.jl)" begin
        doctest(GraphsOptim)
    end

    @testset verbose = true "Min cost flow" begin
        include("minimum_cost_flow.jl")
    end

    @testset verbose = true "Max weight matching" begin
        include("maximum_weight_matching.jl")
    end
end
