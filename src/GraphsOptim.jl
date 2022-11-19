"""
    GraphsOptim

A package for graph-related optimization algorithms that rely on Linear Programming.
"""
module GraphsOptim

using Graphs: AbstractGraph, IsDirected
using Graphs: vertices, edges, nv, ne, src, dst, inneighbors, outneighbors
using FillArrays: Zeros, Fill
using JuMP: Model
using JuMP: objective_function, add_to_expression!
using JuMP: @variable, @constraint, @objective
using SimpleTraits: SimpleTraits, @traitfn

export minimum_cost_flow!

include("minimum_cost_flow.jl")

end
