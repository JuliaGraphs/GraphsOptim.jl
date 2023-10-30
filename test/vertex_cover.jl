using Graphs: SimpleGraph
using GraphsOptim
using Test

adj_matrix = [
    0 1 1 0 0 0
    1 0 1 1 1 1
    1 1 0 0 0 0
    0 1 0 0 0 0
    0 1 0 0 0 0
    0 1 0 0 0 0
]

expected_ans_options = Set([Set([1, 2]), Set([2, 3])])

graph = SimpleGraph(adj_matrix)

@test Set(min_vertex_cover(graph)) in expected_ans_options
