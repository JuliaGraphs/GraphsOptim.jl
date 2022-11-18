"""
    GraphsOptim

A package for graph-related optimization algorithms that rely on Linear Programming.
"""
module GraphsOptim

using Graphs
using JuMP
using SimpleTraits
using SparseArrays

export minimum_cost_flow
export maximum_weight_matching
export maximum_weight_maximal_matching

include("utils.jl")
include("minimum_cost_flow.jl")
include("maximum_weight_matching.jl")
include("maximum_weight_maximal_matching.jl")

end
