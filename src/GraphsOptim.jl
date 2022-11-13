"""
    GraphsOptim

A package for graph-related optimization algorithms that rely on Linear Programming.

The package JuMP.jl and one of its supported solvers are required.
"""
module GraphsOptim

using Graphs
using JuMP
using SimpleTraits
using SparseArrays
using LinearAlgebra

include("utils.jl")
include("minimum_cost_flow.jl")
include("maximum_weight_matching.jl")
include("FAQalgorithm.jl")

export minimum_cost_flow
export maximum_weight_matching, maximum_weight_maximal_matching
export faq, _flat_doubly_stochastic, _gradient, _distance, _solve_transportation_problem,_solve_assignment_problem, _step_size, _update_P

end
