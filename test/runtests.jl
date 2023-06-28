using Aqua
using Documenter
using Graphs
using GraphsOptim
using HiGHS
using JET
using JuliaFormatter
using LinearAlgebra
using SparseArrays
using Test

@testset verbose = true "GraphsOptim.jl" begin
    @testset verbose = true "Quality (Aqua.jl)" begin
        Aqua.test_all(GraphsOptim; ambiguities=false)
    end

    @testset verbose = true "Formatting (JuliaFormatter.jl)" begin
        @test format(GraphsOptim; verbose=false, overwrite=false)
    end

    @testset verbose = true "Correctness (JET.jl)" begin
        if VERSION >= v"1.9"
            JET.test_package(GraphsOptim; target_defined_modules=true)
        end
    end

    @testset verbose = true "Doctests (Documenter.jl)" begin
        doctest(GraphsOptim)
    end

    @testset verbose = true "Minimum cost flow" begin
        include("minimum_cost_flow.jl")
    end

    @testset verbose = true "Alignment" begin
        include("alignment_algorithms.jl")
    end
end;
