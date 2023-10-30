using Graphs
using GraphsOptim
using Test

#=
This graph is the example on the Wikipedia page "Shortest Path Problem"
=#

digraph_adjacency = [
    0 4 2 0 0 0
    0 0 5 10 0 0
    0 0 0 0 3 0
    0 0 0 0 0 11
    0 0 0 4 0 0
    0 0 0 0 0 0
]

undigraph_adjacency = [
    0 4 2 0 0 0
    4 0 5 10 0 0
    2 5 0 0 3 0
    0 10 0 0 4 11
    0 0 3 4 0 0
    0 0 0 11 0 0
]

source, target = 1, 6
answer = [1, 3, 5, 4, 6]

digraph = Graphs.SimpleDiGraph(digraph_adjacency)
digraph_ans = GraphsOptim.shortest_path(digraph, source, target, digraph_adjacency)
@test all(answer .== digraph_ans)

undigraph = Graphs.SimpleGraph(undigraph_adjacency)
undigraph_ans = GraphsOptim.shortest_path(undigraph, source, target, undigraph_adjacency)
@test all(answer .== undigraph_ans)
