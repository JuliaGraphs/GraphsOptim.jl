module GraphsOptim

using Graphs
using JuMP
using SimpleTraits
using SparseArrays

include("mincost_flow.jl")

export mincost_flow

end
