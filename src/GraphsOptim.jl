"""
    GraphsOptim

A package for graph optimization algorithms that rely on mathematical programming.
"""
module GraphsOptim

using Graphs: AbstractGraph, is_directed
using Graphs: vertices, edges, nv, ne, src, dst, inneighbors, outneighbors
using Graphs: complement, maximal_cliques
using FillArrays: Zeros, Fill
using HiGHS: HiGHS
using JuMP: Model
using JuMP: objective_function, add_to_expression!
using JuMP: set_silent, optimize!, termination_status, value
using JuMP: set_optimizer, objective_value
using JuMP: @variable, @constraint, @objective
using LinearAlgebra: norm, tr, dot
using MathOptInterface: OPTIMAL
using SparseArrays: sparse
using OptimalTransport: sinkhorn

export min_cost_flow
export min_cost_assignment
export FAQ, GOAT, graph_matching
export fractional_chromatic_number, fractional_clique_number

include("utils.jl")
include("flow.jl")
include("assignment.jl")
include("graph_matching.jl")
include("fractional_coloring.jl")

end
