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
    @testset "Code formatting" begin
        @test format(GraphsOptim; verbose=false, overwrite=false)
    end

    if VERSION >= v"1.9"
        @testset "Code quality" begin
            Aqua.test_all(GraphsOptim; ambiguities=false)
        end

        @testset "Code linting" begin
            JET.test_package(GraphsOptim; target_defined_modules=true)
        end
    end

    @testset "Doctests" begin
        doctest(GraphsOptim)
    end

    @testset verbose = true "Flow" begin
        include("flow.jl")
    end

    @testset verbose = true "Assignment" begin
        include("assignment.jl")
    end

    @testset verbose = true "Graph matching" begin
        include("graph_matching.jl")
    end

    @testset verbose = true "Vertex cover" begin
        include("min_vertex_cover.jl")
    end

    @testset verbose = true "Cliques" begin
        include("maximum_clique.jl")
    end

    @testset verbose = true "Independent set" begin
        include("independent_set.jl")
    end

    @testset verbose = true "Fractional coloring" begin
        include("fractional_coloring.jl")
    end

    @testset verbose = true "Shortest path" begin
        include("shortest_path.jl")
    end
end;
